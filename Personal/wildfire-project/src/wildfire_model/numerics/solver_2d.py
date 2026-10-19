"""Two-dimensional finite-volume solver for the dry ADR model."""

from dataclasses import dataclass
from math import isfinite
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError, NonPhysicalStateError
from wildfire_model.inputs.fields import (
    FieldInput,
    evaluate_field,
    evaluate_vector_field,
)
from wildfire_model.mesh.mesh_2d import UniformMesh2D
from wildfire_model.models.parameters import ADRParameters
from wildfire_model.models.state import ADRState, validate_temperature_biomass
from wildfire_model.numerics.boundary_conditions import BoundaryCondition2D
from wildfire_model.numerics.cfl import stable_transport_timestep_2d
from wildfire_model.numerics.finite_volume import transport_rhs_2d
from wildfire_model.numerics.reaction_integrator import (
    implicit_midpoint_step,
    implicit_midpoint_step_constant_reaction,
)
from wildfire_model.numerics.time_integration import ssprk3_step_time_dependent
from wildfire_model.numerics.weno import SUPPORTED_WENO_ORDERS

ReactionModel = Literal["arrhenius", "constant"]


@dataclass(frozen=True, slots=True)
class SimulationResult2D:
    """Two-dimensional snapshots with field axes ``(time, y, x)``."""

    times: NDArray[np.float64]
    temperature: NDArray[np.float64]
    biomass: NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class ADRSolver2D:
    """Dry two-dimensional ADR solver based on Reisch Eq. (14).

    Transport uses conservative directional WENO fluxes, central diffusive
    face fluxes, and SSP-RK3. Reaction uses cell-wise implicit midpoint.
    Full steps use the Reisch Sec. 2.2 Strang sequence.
    """

    mesh: UniformMesh2D
    parameters: ADRParameters
    velocity_x: FieldInput = 0.0
    velocity_y: FieldInput = 0.0
    weno_order: int = 7
    boundary: BoundaryCondition2D = "transmissive"
    cfl: float = 0.1
    reaction_model: ReactionModel = "arrhenius"
    constant_reaction_rate: float = 0.01

    def __post_init__(self) -> None:
        if self.weno_order not in SUPPORTED_WENO_ORDERS:
            raise InvalidConfigurationError(f"weno_order must be one of {SUPPORTED_WENO_ORDERS}.")
        if self.boundary not in (
            "periodic",
            "transmissive",
            "zero_diffusive_flux",
            "zero_total_flux",
        ):
            raise InvalidConfigurationError("Unsupported 2D boundary condition.")
        if not isfinite(self.cfl) or not 0.0 < self.cfl <= 1.0:
            raise InvalidConfigurationError("cfl must lie in (0, 1].")
        if self.reaction_model not in ("arrhenius", "constant"):
            raise InvalidConfigurationError("reaction_model must be 'arrhenius' or 'constant'.")
        if not isfinite(self.constant_reaction_rate) or self.constant_reaction_rate < 0.0:
            raise InvalidConfigurationError(
                "constant_reaction_rate must be finite and non-negative."
            )

    def _velocity(self, time: float) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        return evaluate_vector_field(
            self.velocity_x,
            self.velocity_y,
            self.mesh,
            name="transport velocity",
            time=time,
        )

    def maximum_time_step(self, time: float = 0.0) -> float:
        """Return the 2D transport stability limit at ``time`` in s."""

        velocity_x, velocity_y = self._velocity(time)
        return stable_transport_timestep_2d(
            self.mesh,
            velocity_x,
            velocity_y,
            self.parameters.thermal_diffusivity,
            self.cfl,
        )

    def _validate_state_size(self, state: ADRState) -> None:
        if state.temperature.shape != self.mesh.shape:
            raise NonPhysicalStateError(
                f"2D ADR state must have shape {self.mesh.shape}; "
                f"received {state.temperature.shape}."
            )
        validate_temperature_biomass(
            state.temperature,
            state.biomass,
            context="2D ADR solver state",
        )

    def _reaction_step(self, state: ADRState, time_step: float) -> ADRState:
        if self.reaction_model == "arrhenius":
            temperature, biomass = implicit_midpoint_step(
                state.temperature,
                state.biomass,
                time_step,
                self.parameters,
            )
        else:
            temperature, biomass = implicit_midpoint_step_constant_reaction(
                state.temperature,
                state.biomass,
                time_step,
                self.parameters,
                self.constant_reaction_rate,
            )
        return ADRState(temperature, biomass)

    def transport_step(
        self,
        state: ADRState,
        time_step: float,
        *,
        time: float = 0.0,
    ) -> ADRState:
        """Advance two-dimensional transport by one SSP-RK3 step."""

        self._validate_state_size(state)
        if not isfinite(time) or not isfinite(time_step) or time_step <= 0.0:
            raise InvalidConfigurationError(
                "Transport time and time_step must be finite, with time_step > 0 s."
            )
        limits = (
            self.maximum_time_step(time),
            self.maximum_time_step(time + 0.5 * time_step),
            self.maximum_time_step(time + time_step),
        )
        maximum = min(limits)
        if isfinite(maximum) and time_step > maximum * (1.0 + 1.0e-12):
            raise InvalidConfigurationError(
                f"Transport time step {time_step:.16g} s exceeds the 2D CFL limit {maximum:.16g} s."
            )
        diffusivity = evaluate_field(
            self.parameters.thermal_diffusivity,
            self.mesh,
            name="thermal diffusivity",
            minimum=0.0,
        )

        def operator(stage_time: float, values: NDArray[np.float64]) -> NDArray[np.float64]:
            velocity_x, velocity_y = self._velocity(stage_time)
            return transport_rhs_2d(
                values,
                self.mesh,
                velocity_x,
                velocity_y,
                diffusivity,
                self.weno_order,
                self.boundary,
            )

        temperature = ssprk3_step_time_dependent(
            state.temperature,
            time,
            time_step,
            operator,
        )
        validated_temperature, validated_biomass = validate_temperature_biomass(
            temperature,
            state.biomass,
            context="2D transport output",
        )
        return ADRState(validated_temperature, validated_biomass.copy())

    def step(
        self,
        state: ADRState,
        time_step: float,
        *,
        time: float = 0.0,
    ) -> ADRState:
        """Advance one CFL-compliant Strang-split step."""

        first = self._reaction_step(state, 0.5 * time_step)
        transported = self.transport_step(first, time_step, time=time)
        return self._reaction_step(transported, 0.5 * time_step)

    def solve(
        self,
        initial_state: ADRState,
        output_times: ArrayLike,
        *,
        maximum_time_step: float | None = None,
    ) -> SimulationResult2D:
        """Integrate to strictly increasing output times starting at zero."""

        self._validate_state_size(initial_state)
        times = np.asarray(output_times, dtype=np.float64)
        if (
            times.ndim != 1
            or times.size < 1
            or not np.all(np.isfinite(times))
            or times[0] != 0.0
            or np.any(np.diff(times) <= 0.0)
        ):
            raise InvalidConfigurationError(
                "output_times must be finite, increasing, and start at zero."
            )
        if maximum_time_step is not None and (
            not isfinite(maximum_time_step) or maximum_time_step <= 0.0
        ):
            raise InvalidConfigurationError(
                "maximum_time_step must be finite and > 0 s when provided."
            )
        temperature = np.empty((times.size, *self.mesh.shape), dtype=np.float64)
        biomass = np.empty_like(temperature)
        temperature[0] = initial_state.temperature
        biomass[0] = initial_state.biomass
        current = ADRState(
            initial_state.temperature.copy(),
            initial_state.biomass.copy(),
        )
        elapsed = 0.0
        for output_index in range(1, times.size):
            target = float(times[output_index])
            while elapsed < target:
                step_limit = self.maximum_time_step(elapsed)
                if maximum_time_step is not None:
                    step_limit = min(step_limit, maximum_time_step)
                if not isfinite(step_limit):
                    step_limit = target - elapsed
                step = min(step_limit, target - elapsed)
                endpoint_limit = self.maximum_time_step(elapsed + step)
                step = min(step, endpoint_limit)
                current = self.step(current, step, time=elapsed)
                elapsed += step
            temperature[output_index] = current.temperature
            biomass[output_index] = current.biomass
        return SimulationResult2D(times.copy(), temperature, biomass)
