"""Reproducible Navas-Montilla two-dimensional deterministic cases."""

from dataclasses import dataclass, replace
from math import isfinite

import numpy as np
from numpy.typing import NDArray

from wildfire_model.exceptions import InvalidConfigurationError
from wildfire_model.inputs.config import CompleteModelConfiguration
from wildfire_model.inputs.initial_conditions import (
    apply_firebreak,
    circular_ignition,
)
from wildfire_model.mesh.mesh_2d import UniformMesh2D
from wildfire_model.models.state import ADRState


@dataclass(frozen=True, slots=True)
class DeterministicCase2D:
    """Mesh, initial state, physical configuration, and output schedule."""

    mesh: UniformMesh2D
    initial_state: ADRState
    configuration: CompleteModelConfiguration
    output_times: NDArray[np.float64]


def sloping_hill(x: NDArray[np.float64], y: NDArray[np.float64]) -> NDArray[np.float64]:
    """Synthetic terrain in m: a 4% slope and a 6 m Gaussian hill."""

    return 0.04 * x + 6.0 * np.exp(-0.5 * (((x - 26.0) / 9.0) ** 2 + ((y - 16.0) / 10.0) ** 2))


def terrain_fire_case(
    *,
    nx: int = 80,
    ny: int = 64,
    duration: float = 180.0,
    frames: int = 121,
) -> DeterministicCase2D:
    """Build a terrain walkthrough with paper defaults and a synthetic hill."""

    if not isfinite(duration) or duration <= 0.0:
        raise InvalidConfigurationError("duration must be finite and positive in s.")
    if isinstance(frames, bool) or not isinstance(frames, int) or frames < 2:
        raise InvalidConfigurationError("frames must be an integer >= 2.")
    mesh = UniformMesh2D(0.0, 40.0, 0.0, 32.0, nx=nx, ny=ny)
    initial = circular_ignition(
        mesh, (9.0, 16.0), 2.0, ignition_temperature=670.0, ambient_temperature=300.0
    )
    configuration = CompleteModelConfiguration(
        elevation=sloping_hill,
        wind_x=1.0,
        moisture_model="S2M",
        moisture=0.1,
    )
    return DeterministicCase2D(mesh, initial, configuration, np.linspace(0.0, duration, frames))


def _paper_mesh(full_resolution: bool) -> UniformMesh2D:
    cells = 400 if full_resolution else 80
    return UniformMesh2D(0.0, 200.0, 0.0, 200.0, cells, cells)


def homogeneous_moisture_case(
    moisture: float,
    *,
    full_resolution: bool = False,
) -> DeterministicCase2D:
    """Return Navas-Montilla Sec. 4.5.1, Eqs. (80)-(81)."""

    mesh = _paper_mesh(full_resolution)
    state = circular_ignition(
        mesh,
        (100.0, 100.0),
        5.0,
        ignition_temperature=670.0,
        ambient_temperature=300.0,
    )
    configuration = CompleteModelConfiguration(
        moisture_model="S2M",
        moisture=moisture,
        wind_x=0.25,
        wind_y=0.25,
        combustion_model="constant",
        mixture_approximation="fuel_only",
    )
    final_time = 1400.0 if full_resolution else 200.0
    return DeterministicCase2D(
        mesh,
        state,
        configuration,
        np.linspace(0.0, final_time, 8 if full_resolution else 5),
    )


def quadrant_moisture_case(
    *,
    full_resolution: bool = False,
) -> DeterministicCase2D:
    """Return Navas-Montilla Sec. 4.5.2, Eqs. (82)-(83)."""

    mesh = _paper_mesh(full_resolution)
    x, y = mesh.cell_centers
    moisture = np.where(
        y > 100.0,
        np.where(x < 100.0, 0.0, 0.1),
        np.where(x < 100.0, 0.2, 0.3),
    )
    state = circular_ignition(
        mesh,
        (100.0, 100.0),
        5.0,
        ignition_temperature=670.0,
        ambient_temperature=300.0,
    )
    configuration = CompleteModelConfiguration(
        moisture_model="S2M",
        moisture=moisture,
        combustion_model="constant",
        mixture_approximation="fuel_only",
    )
    final_time = 1400.0 if full_resolution else 200.0
    return DeterministicCase2D(
        mesh,
        state,
        configuration,
        np.linspace(0.0, final_time, 8 if full_resolution else 5),
    )


def heterogeneous_coverage_case(
    reduced_coverage: float,
    *,
    full_resolution: bool = False,
) -> DeterministicCase2D:
    """Return Navas-Montilla Sec. 4.5.3, Eqs. (84)-(85)."""

    mesh = _paper_mesh(full_resolution)
    x, _ = mesh.cell_centers
    coverage = np.where(x <= 100.0, 1.0, reduced_coverage)
    state = circular_ignition(
        mesh,
        (50.0, 50.0),
        5.0,
        ignition_temperature=670.0,
        ambient_temperature=300.0,
    )
    configuration = CompleteModelConfiguration(
        surface_coverage=coverage,
        moisture_model="S2M",
        moisture=0.1,
        wind_x=1.0,
        wind_y=1.0,
        combustion_model="constant",
        mixture_approximation="fuel_only",
    )
    final_time = 2200.0 if full_resolution else 250.0
    return DeterministicCase2D(
        mesh,
        state,
        configuration,
        np.linspace(0.0, final_time, 12 if full_resolution else 6),
    )


def topography_case(
    *,
    include_topography: bool,
    full_resolution: bool = False,
) -> DeterministicCase2D:
    """Return Navas-Montilla Sec. 4.5.4, Eqs. (86)-(88)."""

    mesh = _paper_mesh(full_resolution)
    x, y = mesh.cell_centers
    fuel_load = (
        np.sin(x / 20.0) * np.cos(y / 10.0 - x / 20.0)
        + np.sin(x / 20.0 - 2.0) * np.sin(y / 10.0)
        + 2.0
    )
    elevation = (
        50.0 * np.exp(-((x - 200.0) ** 2 + (y - 200.0) ** 2) / 20_000.0)
        if include_topography
        else np.zeros(mesh.shape)
    )
    state = circular_ignition(
        mesh,
        (100.0, 50.0),
        5.0,
        ignition_temperature=670.0,
        ambient_temperature=300.0,
    )
    state = ADRState.from_arrays(
        state.temperature,
        apply_firebreak(
            state.biomass,
            mesh,
            fuel_load <= 0.0,
            remaining_fraction=0.0,
        ),
    )
    configuration = CompleteModelConfiguration(
        maximum_fuel_load=fuel_load,
        moisture_model="NM",
        wind_y=0.5,
        elevation=elevation,
        combustion_model="constant",
        mixture_approximation="complete",
    )
    final_time = 2600.0 if full_resolution else 300.0
    return DeterministicCase2D(
        mesh,
        state,
        configuration,
        np.linspace(0.0, final_time, 14 if full_resolution else 7),
    )


def with_arrhenius_combustion(
    case: DeterministicCase2D,
) -> DeterministicCase2D:
    """Return a case using the primary Arrhenius combustion option."""

    return replace(
        case,
        configuration=replace(
            case.configuration,
            combustion_model="arrhenius",
        ),
    )
