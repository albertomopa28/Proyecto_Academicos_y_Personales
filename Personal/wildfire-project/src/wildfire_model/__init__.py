"""Deterministic advection-diffusion-reaction wildfire model."""

from wildfire_model.helpers.paraview import write_paraview_series
from wildfire_model.inputs.config import CompleteModelConfiguration
from wildfire_model.mesh.mesh_1d import UniformMesh1D
from wildfire_model.mesh.mesh_2d import UniformMesh2D
from wildfire_model.models.parameters import ADRParameters
from wildfire_model.models.state import ADRState
from wildfire_model.numerics.complete_solver_2d import NavasMontillaSolver2D
from wildfire_model.numerics.solver_1d import ADRSolver1D, SimulationResult
from wildfire_model.numerics.solver_2d import ADRSolver2D, SimulationResult2D

__all__ = [
    "ADRParameters",
    "ADRSolver1D",
    "ADRSolver2D",
    "ADRState",
    "CompleteModelConfiguration",
    "NavasMontillaSolver2D",
    "SimulationResult",
    "SimulationResult2D",
    "UniformMesh1D",
    "UniformMesh2D",
    "write_paraview_series",
]
