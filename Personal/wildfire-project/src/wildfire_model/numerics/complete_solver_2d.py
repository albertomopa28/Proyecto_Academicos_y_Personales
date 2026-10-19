"""Finite-volume solver for the deterministic Navas-Montilla model."""

from dataclasses import dataclass, field
from math import isfinite

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError, NonPhysicalStateError
from wildfire_model.inputs.config import CompleteModelConfiguration
from wildfire_model.inputs.fields import evaluate_field, evaluate_vector_field
from wildfire_model.mesh.mesh_2d import UniformMesh2D
from wildfire_model.models.moisture import apparent_specific_heat
from wildfire_model.models.properties import (
    effective_thermal_conductivity,
    effective_velocity,
    initial_fuel_load,
    initial_fuel_volume_fraction,
    mixture_properties,
)
from wildfire_model.models.state import ADRState, validate_temperature_biomass
from wildfire_model.numerics.cfl import stable_transport_timestep_2d
from wildfire_model.numerics.finite_volume import energy_transport_rhs_2d
from wildfire_model.numerics.solver_2d import SimulationResult2D
from wildfire_model.numerics.time_integration import ssprk3_step_time_dependent
from wildfire_model.numerics.weno import SUPPORTED_WENO_ORDERS


@dataclass(frozen=True, slots=True)
class NavasMontillaSolver2D:
    """Solve Navas-Montilla Sec. 2.2, Eq. (32), with SSP-RK3.

    The apparent-capacity source and nonlinear conductivity are evaluated at
    every Runge-Kutta stage. C2M remains unavailable because the paper does
    not define the post-evaporation heat capacity or Eq. (67) SI units.
    """

    mesh: UniformMesh2D
    configuration: CompleteModelConfiguration = field(default_factory=CompleteModelConfiguration)

    def __post_init__(self) -> None:
        if self.configuration.weno_order not in SUPPORTED_WENO_ORDERS:
            raise InvalidConfigurationError(f"weno_order must be one of {SUPPORTED_WENO_ORDERS}.")
        if not isfinite(self.configuration.cfl) or not 0.0 < self.configuration.cfl <= 1.0:
            raise InvalidConfigurationError("cfl must lie in (0, 1].")
        if (
            not isfinite(self.configuration.arrhenius_pre_exponential_factor)
            or self.configuration.arrhenius_pre_exponential_factor < 0.0
            or not isfinite(self.configuration.constant_reaction_rate)
            or self.configuration.constant_reaction_rate < 0.0
        ):
            raise InvalidConfigurationError("Combustion rates must be finite and non-negative.")

    def _spatial_inputs(
        self,
        time: float,
    ) -> tuple[
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
    ]:
        config = self.configuration
        coverage = evaluate_field(
            config.surface_coverage,
            self.mesh,
            name="surface coverage",
            time=time,
            minimum=0.0,
            maximum=1.0,
        )
        maximum_load = evaluate_field(
            (
                config.fuel.maximum_fuel_load
                if config.maximum_fuel_load is None
                else config.maximum_fuel_load
            ),
            self.mesh,
            name="maximum fuel load",
            time=time,
            minimum=0.0,
        )
        moisture = evaluate_field(
            config.moisture,
            self.mesh,
            name="fuel moisture",
            time=time,
            minimum=0.0,
        )
        conductivity = evaluate_field(
            (
                config.calibration.conductive_thermal_conductivity
                if config.conductive_thermal_conductivity is None
                else config.conductive_thermal_conductivity
            ),
            self.mesh,
            name="conductive thermal conductivity",
            time=time,
            minimum=0.0,
        )
        optical_path = evaluate_field(
            (
                config.calibration.optical_path_length
                if config.optical_path_length is None
                else config.optical_path_length
            ),
            self.mesh,
            name="optical path length",
            time=time,
            minimum=0.0,
        )
        ambient = evaluate_field(
            (
                config.environment.ambient_temperature
                if config.ambient_temperature is None
                else config.ambient_temperature
            ),
            self.mesh,
            name="ambient temperature",
            time=time,
            minimum=np.finfo(np.float64).tiny,
        )
        wind_x, wind_y = evaluate_vector_field(
            config.wind_x,
            config.wind_y,
            self.mesh,
            name="wind",
            time=time,
        )
        elevation = evaluate_field(
            config.elevation,
            self.mesh,
            name="elevation",
            time=time,
        )
        velocity_x, velocity_y = effective_velocity(
            wind_x,
            wind_y,
            elevation,
            self.mesh,
            wind_correction=config.calibration.wind_correction,
            slope_correction=config.calibration.slope_correction,
        )
        return (
            coverage,
            maximum_load,
            moisture,
            conductivity,
            optical_path,
            ambient,
            velocity_x,
            velocity_y,
        )

    def _physical_fields(
        self,
        temperature: NDArray[np.float64],
        biomass: NDArray[np.float64],
        time: float,
    ) -> tuple[
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
        NDArray[np.float64],
    ]:
        config = self.configuration
        (
            coverage,
            maximum_load,
            moisture,
            conductive,
            optical_path,
            ambient,
            velocity_x,
            velocity_y,
        ) = self._spatial_inputs(time)
        load = initial_fuel_load(maximum_load, coverage)
        initial_fraction = initial_fuel_volume_fraction(
            load,
            config.fuel.intrinsic_density,
            config.fuel.layer_height,
        )
        fuel_specific_heat = apparent_specific_heat(
            config.moisture_model,
            temperature,
            biomass,
            moisture,
            dry_specific_heat=config.fuel.dry_specific_heat,
            water_specific_heat=config.environment.water_specific_heat,
            latent_heat=config.environment.water_latent_heat,
            ambient_temperature=config.environment.ambient_temperature,
            evaporation_temperature=config.environment.water_evaporation_temperature,
            pyrolysis_temperature=config.fuel.pyrolysis_temperature,
            evaporation_interval=config.evaporation_interval,
        )
        density, mixture_heat = mixture_properties(
            config.fuel.intrinsic_density,
            fuel_specific_heat,
            config.environment.air_density,
            config.environment.air_specific_heat,
            initial_fraction,
            biomass,
            config.mixture_approximation,
        )
        capacity = density * mixture_heat
        thermal_conductivity = effective_thermal_conductivity(
            temperature,
            conductive,
            config.environment.emissivity,
            optical_path,
            config.environment.stefan_boltzmann_constant,
        )
        return (
            initial_fraction,
            capacity,
            thermal_conductivity,
            ambient,
            velocity_x,
            velocity_y,
            moisture,
        )

    def _combustion_rate(self, temperature: NDArray[np.float64]) -> NDArray[np.float64]:
        config = self.configuration
        active = temperature >= config.fuel.pyrolysis_temperature
        if config.combustion_model == "arrhenius":
            return np.where(
                active,
                config.arrhenius_pre_exponential_factor
                * np.exp(-config.fuel.pyrolysis_temperature / temperature),
                0.0,
            )
        return np.where(active, config.constant_reaction_rate, 0.0)

    def right_hand_side(
        self,
        time: float,
        state: NDArray[np.float64],
    ) -> NDArray[np.float64]:
        """Return the pure semi-discrete residual for ``(T, Y)``."""

        if state.shape != (2, *self.mesh.shape):
            raise InvalidConfigurationError(
                f"Complete-model state must have shape {(2, *self.mesh.shape)}."
            )
        temperature, biomass = validate_temperature_biomass(
            state[0],
            state[1],
            context="complete 2D RK stage",
        )
        (
            initial_fraction,
            capacity,
            conductivity,
            ambient,
            velocity_x,
            velocity_y,
            _,
        ) = self._physical_fields(temperature, biomass, time)
        transport = energy_transport_rhs_2d(
            temperature,
            self.mesh,
            velocity_x,
            velocity_y,
            conductivity,
            capacity,
            self.configuration.weno_order,
            self.configuration.boundary,
        )
        rate = self._combustion_rate(temperature)
        biomass_rate = -rate * biomass
        cooling_coefficient = (
            self.configuration.calibration.surface_cooling_coefficient
            / self.configuration.fuel.layer_height
        )
        heat_source = (
            self.configuration.fuel.intrinsic_density
            * initial_fraction
            * self.configuration.fuel.heating_value
            * rate
            * biomass
            / capacity
        )
        temperature_rate = (
            transport - cooling_coefficient * (temperature - ambient) / capacity + heat_source
        )
        return np.stack((temperature_rate, biomass_rate))

    def maximum_time_step(self, state: ADRState, time: float = 0.0) -> float:
        """Return the combined transport and explicit-reaction step limit."""

        self._validate_state(state)
        (
            _,
            capacity,
            conductivity,
            _,
            velocity_x,
            velocity_y,
            _,
        ) = self._physical_fields(state.temperature, state.biomass, time)
        transport_limit = stable_transport_timestep_2d(
            self.mesh,
            velocity_x,
            velocity_y,
            conductivity / capacity,
            self.configuration.cfl,
        )
        maximum_rate = float(np.max(self._combustion_rate(state.temperature)))
        reaction_limit = (
            float("inf") if maximum_rate == 0.0 else self.configuration.cfl / maximum_rate
        )
        return min(transport_limit, reaction_limit)

    def _validate_state(self, state: ADRState) -> None:
        if state.temperature.shape != self.mesh.shape:
            raise NonPhysicalStateError(f"Complete-model fields must have shape {self.mesh.shape}.")
        validate_temperature_biomass(
            state.temperature,
            state.biomass,
            context="complete 2D model",
        )

    def step(
        self,
        state: ADRState,
        time_step: float,
        *,
        time: float = 0.0,
    ) -> ADRState:
        """Advance Eq. (32) with Navas-Montilla Eq. (78)."""

        self._validate_state(state)
        if not isfinite(time_step) or time_step <= 0.0:
            raise InvalidConfigurationError("time_step must be finite and positive.")
        maximum = min(
            self.maximum_time_step(state, time),
            self.maximum_time_step(state, time + time_step),
        )
        if isfinite(maximum) and time_step > maximum * (1.0 + 1.0e-12):
            raise InvalidConfigurationError(
                f"Time step {time_step:.16g} s exceeds the complete-model limit {maximum:.16g} s."
            )
        packed = np.stack((state.temperature, state.biomass))
        result = ssprk3_step_time_dependent(
            packed,
            time,
            time_step,
            self.right_hand_side,
        )
        temperature, biomass = validate_temperature_biomass(
            result[0],
            result[1],
            context="complete 2D step output",
        )
        return ADRState(temperature, biomass)

    def solve(
        self,
        initial_state: ADRState,
        output_times: ArrayLike,
        *,
        maximum_time_step: float | None = None,
    ) -> SimulationResult2D:
        """Integrate the complete model to each specified output time."""

        self._validate_state(initial_state)
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
            raise InvalidConfigurationError("maximum_time_step must be finite and positive.")
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
                limit = self.maximum_time_step(current, elapsed)
                if maximum_time_step is not None:
                    limit = min(limit, maximum_time_step)
                if not isfinite(limit):
                    limit = target - elapsed
                step = min(limit, target - elapsed)
                while True:
                    try:
                        current = self.step(current, step, time=elapsed)
                    except NonPhysicalStateError:
                        step *= 0.5
                        if step <= np.finfo(np.float64).eps * max(1.0, abs(target)):
                            raise
                    else:
                        break
                elapsed += step
            temperature[output_index] = current.temperature
            biomass[output_index] = current.biomass
        return SimulationResult2D(times.copy(), temperature, biomass)
