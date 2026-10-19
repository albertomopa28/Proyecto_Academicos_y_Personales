"""Finite-volume numerical operators."""

from wildfire_model.numerics.solver_1d import ADRSolver1D, SimulationResult
from wildfire_model.numerics.solver_2d import ADRSolver2D, SimulationResult2D

__all__ = [
    "ADRSolver1D",
    "ADRSolver2D",
    "SimulationResult",
    "SimulationResult2D",
]
