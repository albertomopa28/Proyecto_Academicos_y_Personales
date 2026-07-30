"""One-dimensional finite-volume composition of the base ADR model."""

from dataclasses import dataclass
from math import isfinite

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError, NonPhysicalStateError
from wildfire_model.mesh.mesh_1d import UniformMesh1D
from wildfire_model.models.parameters import ADRParameters
from wildfire_model.models.state import ADRState, validate_temperature_biomass
from wildfire_model.numerics.boundary_conditions import BoundaryCondition
from wildfire_model.numerics.cfl import stable_transport_timestep
from wildfire_model.numerics.finite_volume import transport_rhs
from wildfire_model.numerics.reaction_integrator import implicit_midpoint_step
from wildfire_model.numerics.time_integration import ssprk3_step
from wildfire_model.numerics.weno import SUPPORTED_WENO_ORDERS


@dataclass(frozen=True, slots=True)
class SimulationResult:
    """Snapshots returned by :meth:`ADRSolver1D.solve`."""

    times: NDArray[np.float64]
    temperature: NDArray[np.float64]
    biomass: NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class ADRSolver1D:
    """Stage-1 solver for Reisch et al. (2024), Eq. (14).

    Advection-diffusion is advanced by SSP-RK3. The local reaction and
    ambient cooling system is advanced by implicit midpoint. A full step
    is reaction(dt/2), transport(dt), reaction(dt/2), i.e. Strang
    splitting as specified in Reisch et al. Sec. 2.2.
    """

    mesh: UniformMesh1D
    parameters: ADRParameters
    weno_order: int = 7
    boundary: BoundaryCondition = "transmissive"
    cfl: float = 0.1

    def __post_init__(self) -> None:
        if self.weno_order not in SUPPORTED_WENO_ORDERS:
            raise InvalidConfigurationError(f"weno_order must be one of {SUPPORTED_WENO_ORDERS}.")
        if self.boundary not in ("periodic", "transmissive"):
            raise InvalidConfigurationError("boundary must be either 'periodic' or 'transmissive'.")
        if not isfinite(self.cfl) or not 0.0 < self.cfl <= 1.0:
            raise InvalidConfigurationError("cfl must lie in (0, 1].")

    @property
    def maximum_time_step(self) -> float:
        """Transport CFL limit in s."""

        return stable_transport_timestep(
            self.mesh,
            self.parameters.velocity,
            self.parameters.thermal_diffusivity,
            self.cfl,
        )

    def _validate_state_size(self, state: ADRState) -> None:
        if state.temperature.shape != (self.mesh.n_cells,):
            raise NonPhysicalStateError(
                "ADR state must contain exactly one cell average per mesh cell; "
                f"expected {self.mesh.n_cells}, received {state.temperature.shape}."
            )
        validate_temperature_biomass(state.temperature, state.biomass, context="ADR solver state")

    def transport_step(self, state: ADRState, time_step: float) -> ADRState:
        """Advance only advection and diffusion for one full step."""

        self._validate_state_size(state)
        if not isfinite(time_step) or time_step <= 0.0:
            raise InvalidConfigurationError("time_step must be finite and > 0 s.")
        maximum = self.maximum_time_step
        if isfinite(maximum) and time_step > maximum * (1.0 + 1.0e-12):
            raise InvalidConfigurationError(
                f"Transport time step {time_step:.16g} s exceeds the CFL limit {maximum:.16g} s."
            )

        def operator(values: NDArray[np.float64]) -> NDArray[np.float64]:
            return transport_rhs(
                values,
                self.mesh,
                self.parameters.velocity,
                self.parameters.thermal_diffusivity,
                self.weno_order,
                self.boundary,
            )

        temperature = ssprk3_step(state.temperature, time_step, operator)
        validated_temperature, validated_biomass = validate_temperature_biomass(
            temperature, state.biomass, context="transport output"
        )
        return ADRState(validated_temperature, validated_biomass.copy())

    def step(self, state: ADRState, time_step: float) -> ADRState:
        """Advance one CFL-compliant Strang-split step."""

        self._validate_state_size(state)
        if not isfinite(time_step) or time_step <= 0.0:
            raise InvalidConfigurationError("time_step must be finite and > 0 s.")
        maximum = self.maximum_time_step
        if isfinite(maximum) and time_step > maximum * (1.0 + 1.0e-12):
            raise InvalidConfigurationError(
                f"Strang time step {time_step:.16g} s exceeds the CFL limit {maximum:.16g} s."
            )
        first_temperature, first_biomass = implicit_midpoint_step(
            state.temperature,
            state.biomass,
            0.5 * time_step,
            self.parameters,
        )
        transported = self.transport_step(ADRState(first_temperature, first_biomass), time_step)
        final_temperature, final_biomass = implicit_midpoint_step(
            transported.temperature,
            transported.biomass,
            0.5 * time_step,
            self.parameters,
        )
        return ADRState(final_temperature, final_biomass)

    def solve(
        self,
        initial_state: ADRState,
        output_times: ArrayLike,
        *,
        maximum_time_step: float | None = None,
    ) -> SimulationResult:
        """Integrate to each strictly increasing output time.

        The first output time must be zero. Steps are shortened exactly
        at output times; otherwise they use the lesser of the CFL limit
        and the optional configured limit.
        """

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
                "output_times must be a finite, strictly increasing 1D array starting at 0."
            )
        step_limit = self.maximum_time_step
        if maximum_time_step is not None:
            if not isfinite(maximum_time_step) or maximum_time_step <= 0.0:
                raise InvalidConfigurationError(
                    "maximum_time_step must be finite and > 0 s when provided."
                )
            step_limit = min(step_limit, maximum_time_step)
        if not isfinite(step_limit):
            step_limit = float(times[-1]) if times[-1] > 0.0 else 1.0

        temperatures = np.empty((times.size, self.mesh.n_cells), dtype=np.float64)
        biomass = np.empty_like(temperatures)
        temperatures[0] = initial_state.temperature
        biomass[0] = initial_state.biomass
        current = ADRState(
            initial_state.temperature.copy(),
            initial_state.biomass.copy(),
        )
        elapsed = 0.0
        for output_index in range(1, times.size):
            target = float(times[output_index])
            while elapsed < target:
                time_step = min(step_limit, target - elapsed)
                current = self.step(current, time_step)
                elapsed += time_step
            temperatures[output_index] = current.temperature
            biomass[output_index] = current.biomass
        return SimulationResult(times.copy(), temperatures, biomass)
