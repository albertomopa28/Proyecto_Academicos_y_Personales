"""TOML configuration loading for one-dimensional ADR simulations."""

import tomllib
from dataclasses import dataclass, field
from math import isfinite
from pathlib import Path
from typing import cast

import numpy as np
from numpy.typing import NDArray

from wildfire_model.exceptions import InvalidConfigurationError
from wildfire_model.inputs.fields import FieldInput
from wildfire_model.models.combustion import CombustionModel
from wildfire_model.models.moisture import MoistureModel
from wildfire_model.models.parameters import (
    ADRParameters,
    CalibrationParameters,
    EnvironmentalProperties,
    FuelProperties,
)
from wildfire_model.models.properties import MixtureApproximation
from wildfire_model.numerics.boundary_conditions import (
    BoundaryCondition,
    BoundaryCondition2D,
)


@dataclass(frozen=True, slots=True)
class MeshConfiguration:
    """Uniform one-dimensional mesh input."""

    x_left: float
    x_right: float
    n_cells: int


@dataclass(frozen=True, slots=True)
class SimulationConfiguration:
    """Complete deterministic one-dimensional simulation input."""

    parameters: ADRParameters
    mesh: MeshConfiguration
    cfl: float
    weno_order: int
    boundary: BoundaryCondition
    output_times: NDArray[np.float64]
    maximum_time_step: float | None
    initial_condition: str


@dataclass(frozen=True, slots=True)
class CompleteModelConfiguration:
    """Separated physical and spatial inputs for Navas-Montilla Eq. (32)."""

    fuel: FuelProperties = field(default_factory=FuelProperties)
    environment: EnvironmentalProperties = field(default_factory=EnvironmentalProperties)
    calibration: CalibrationParameters = field(default_factory=CalibrationParameters)
    surface_coverage: FieldInput = 1.0
    maximum_fuel_load: FieldInput | None = None
    moisture: FieldInput = 0.0
    wind_x: FieldInput = 0.0
    wind_y: FieldInput = 0.0
    elevation: FieldInput = 0.0
    conductive_thermal_conductivity: FieldInput | None = None
    optical_path_length: FieldInput | None = None
    ambient_temperature: FieldInput | None = None
    mixture_approximation: MixtureApproximation = "complete"
    moisture_model: MoistureModel = "NM"
    combustion_model: CombustionModel = "arrhenius"
    arrhenius_pre_exponential_factor: float = 0.0173
    constant_reaction_rate: float = 0.01
    evaporation_interval: float | None = None
    boundary: BoundaryCondition2D = "transmissive"
    weno_order: int = 7
    cfl: float = 0.4


def load_configuration(path: str | Path) -> SimulationConfiguration:
    """Load physical, numerical, mesh, and simulation inputs from TOML."""

    with Path(path).open("rb") as stream:
        data = tomllib.load(stream)
    physics = data["physics"]
    numerics = data["numerics"]
    mesh = data["mesh"]
    simulation = data["simulation"]

    parameters = ADRParameters(
        bulk_density=float(physics["bulk_density"]),
        specific_heat=float(physics["specific_heat"]),
        pre_exponential_factor=float(physics["pre_exponential_factor"]),
        heating_value=float(physics["heating_value"]),
        cooling_coefficient=float(physics["cooling_coefficient"]),
        ambient_temperature=float(physics["ambient_temperature"]),
        pyrolysis_temperature=float(physics["pyrolysis_temperature"]),
        activation_temperature=float(physics["activation_temperature"]),
        thermal_conductivity=float(physics["thermal_conductivity"]),
        velocity=float(physics["velocity"]),
    )
    final_time = float(simulation["final_time"])
    output_interval = float(simulation["output_interval"])
    if not isfinite(final_time) or final_time < 0.0:
        raise InvalidConfigurationError("final_time must be finite and non-negative.")
    if not isfinite(output_interval) or output_interval <= 0.0:
        raise InvalidConfigurationError("output_interval must be finite and positive.")
    output_times = np.arange(0.0, final_time, output_interval, dtype=np.float64)
    output_times = np.append(output_times, final_time)
    maximum_time_step_value = simulation.get("maximum_time_step")
    maximum_time_step = None if maximum_time_step_value is None else float(maximum_time_step_value)
    boundary_value = str(numerics["boundary"])
    if boundary_value not in ("periodic", "transmissive"):
        raise ValueError("boundary must be either 'periodic' or 'transmissive'.")
    boundary = cast(BoundaryCondition, boundary_value)
    return SimulationConfiguration(
        parameters=parameters,
        mesh=MeshConfiguration(
            x_left=float(mesh["x_left"]),
            x_right=float(mesh["x_right"]),
            n_cells=int(mesh["n_cells"]),
        ),
        cfl=float(numerics["cfl"]),
        weno_order=int(numerics["weno_order"]),
        boundary=boundary,
        output_times=np.asarray(output_times, dtype=np.float64),
        maximum_time_step=maximum_time_step,
        initial_condition=str(simulation["initial_condition"]),
    )
