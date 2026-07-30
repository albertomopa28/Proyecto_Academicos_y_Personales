"""Public entry point for deterministic one-dimensional ADR simulations."""

from dataclasses import dataclass
from pathlib import Path

from wildfire_model.exceptions import InvalidConfigurationError
from wildfire_model.inputs.config import SimulationConfiguration, load_configuration
from wildfire_model.inputs.reference_cases import reisch_wave_initial_state
from wildfire_model.mesh.mesh_1d import UniformMesh1D
from wildfire_model.models.parameters import ADRParameters
from wildfire_model.numerics.solver_1d import ADRSolver1D, SimulationResult


@dataclass(frozen=True, slots=True)
class SimulationRun:
    """Mesh, physical parameters, and snapshots from a simulation."""

    mesh: UniformMesh1D
    parameters: ADRParameters
    result: SimulationResult


def run_simulation(configuration: SimulationConfiguration) -> SimulationRun:
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
    return SimulationRun(mesh, configuration.parameters, result)


def run_from_config(path: str | Path) -> SimulationRun:
    """Load a TOML configuration and run its simulation."""

    return run_simulation(load_configuration(path))
