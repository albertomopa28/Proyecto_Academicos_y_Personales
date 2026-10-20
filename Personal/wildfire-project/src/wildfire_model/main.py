"""Simulation entry points with ParaView output under downloads/."""

from dataclasses import dataclass
from pathlib import Path

from wildfire_model.exceptions import InvalidConfigurationError
from wildfire_model.helpers.paraview import write_paraview_series
from wildfire_model.inputs.config import SimulationConfiguration, load_configuration
from wildfire_model.inputs.reference_cases import reisch_wave_initial_state
from wildfire_model.inputs.reference_cases_2d import DeterministicCase2D
from wildfire_model.mesh.mesh_1d import UniformMesh1D
from wildfire_model.models.parameters import ADRParameters
from wildfire_model.numerics.complete_solver_2d import NavasMontillaSolver2D
from wildfire_model.numerics.solver_1d import ADRSolver1D, SimulationResult
from wildfire_model.numerics.solver_2d import SimulationResult2D


@dataclass(frozen=True, slots=True)
class SimulationRun:
    """Mesh, physical parameters, and snapshots from a simulation."""

    mesh: UniformMesh1D
    parameters: ADRParameters
    result: SimulationResult


def run_simulation(
    configuration: SimulationConfiguration,
    *,
    output_directory: str | Path | None = "downloads/reisch_1d",
) -> SimulationRun:
    """Create the mesh, model inputs, solver, and execute the simulation."""

    mesh = UniformMesh1D(
        configuration.mesh.x_left,
        configuration.mesh.x_right,
        configuration.mesh.n_cells,
    )
    if configuration.initial_condition != "reisch_wave":
        raise InvalidConfigurationError(
            f"Unsupported initial condition {configuration.initial_condition!r}."
        )
    initial_state = reisch_wave_initial_state(mesh)
    solver = ADRSolver1D(
        mesh,
        configuration.parameters,
        configuration.weno_order,
        configuration.boundary,
        configuration.cfl,
    )
    result = solver.solve(
        initial_state,
        configuration.output_times,
        maximum_time_step=configuration.maximum_time_step,
    )
    if output_directory is not None:
        write_paraview_series(mesh, result, output_directory=output_directory)
    return SimulationRun(mesh, configuration.parameters, result)


def run_from_config(
    path: str | Path,
    *,
    output_directory: str | Path | None = "downloads/reisch_1d",
) -> SimulationRun:
    """Load a TOML configuration and run its simulation."""

    return run_simulation(load_configuration(path), output_directory=output_directory)


def run_complete_case(
    case: DeterministicCase2D,
    *,
    maximum_time_step: float | None = None,
    output_directory: str | Path | None = "downloads/wildfire_2d",
) -> SimulationResult2D:
    """Run a reproducible two-dimensional Navas-Montilla case."""

    solver = NavasMontillaSolver2D(case.mesh, case.configuration)
    result = solver.solve(
        case.initial_state,
        case.output_times,
        maximum_time_step=maximum_time_step,
    )
    if output_directory is not None:
        write_paraview_series(
            case.mesh,
            result,
            output_directory=output_directory,
            elevation=case.configuration.elevation,
        )
    return result
