"""Uniform one-dimensional finite-volume mesh."""

from dataclasses import dataclass
from math import isfinite

import numpy as np
from numpy.typing import NDArray

from wildfire_model.exceptions import InvalidConfigurationError


@dataclass(frozen=True, slots=True)
class UniformMesh1D:
    """Cell-centred Cartesian mesh on ``[x_left, x_right]``.

    Cell ``i`` occupies ``[x_left + i*dx, x_left + (i+1)*dx]``.
    Stored solution values are cell averages, and ``cell_centers[i]`` is
    the geometric centre of that cell.
    """

    x_left: float
    x_right: float
    n_cells: int

    def __post_init__(self) -> None:
        if not isfinite(self.x_left) or not isfinite(self.x_right):
            raise InvalidConfigurationError("Mesh bounds must be finite.")
        if self.x_right <= self.x_left:
            raise InvalidConfigurationError("Mesh requires x_right > x_left.")
        if isinstance(self.n_cells, bool) or self.n_cells < 1:
            raise InvalidConfigurationError("Mesh n_cells must be a positive integer.")

    @property
    def length(self) -> float:
        """Domain length in m."""

        return self.x_right - self.x_left

    @property
    def spacing(self) -> float:
        """Uniform cell width ``dx`` in m."""

        return self.length / self.n_cells

    @property
    def cell_centers(self) -> NDArray[np.float64]:
        """Cell-centre coordinates in m."""

        return self.x_left + (np.arange(self.n_cells, dtype=np.float64) + 0.5) * self.spacing

    @property
    def faces(self) -> NDArray[np.float64]:
        """Cell-face coordinates in m."""

        return np.linspace(self.x_left, self.x_right, self.n_cells + 1, dtype=np.float64)
