"""Uniform two-dimensional Cartesian finite-volume mesh."""

from dataclasses import dataclass
from math import isfinite

import numpy as np
from numpy.typing import NDArray

from wildfire_model.exceptions import InvalidConfigurationError


@dataclass(frozen=True, slots=True)
class UniformMesh2D:
    """Cell-centred mesh with fields stored as ``(ny, nx)`` arrays.

    Axis 0 is the y coordinate and axis 1 is the x coordinate. Cell
    ``(j, i)`` occupies ``[x_min+i*dx, x_min+(i+1)*dx]`` by
    ``[y_min+j*dy, y_min+(j+1)*dy]``.
    """

    x_min: float
    x_max: float
    y_min: float
    y_max: float
    nx: int
    ny: int

    def __post_init__(self) -> None:
        bounds = (self.x_min, self.x_max, self.y_min, self.y_max)
        if not all(isfinite(value) for value in bounds):
            raise InvalidConfigurationError("Mesh bounds must be finite.")
        if self.x_max <= self.x_min or self.y_max <= self.y_min:
            raise InvalidConfigurationError("Mesh requires x_max > x_min and y_max > y_min.")
        if any(
            isinstance(count, bool) or not isinstance(count, (int, np.integer)) or count < 1
            for count in (self.nx, self.ny)
        ):
            raise InvalidConfigurationError("Mesh nx and ny must be positive integers.")

    @property
    def shape(self) -> tuple[int, int]:
        """Canonical field shape ``(ny, nx)``."""

        return self.ny, self.nx

    @property
    def length_x(self) -> float:
        """Domain length in the x direction in m."""

        return self.x_max - self.x_min

    @property
    def length_y(self) -> float:
        """Domain length in the y direction in m."""

        return self.y_max - self.y_min

    @property
    def dx(self) -> float:
        """Uniform x-directed cell width in m."""

        return self.length_x / self.nx

    @property
    def dy(self) -> float:
        """Uniform y-directed cell width in m."""

        return self.length_y / self.ny

    @property
    def cell_area(self) -> float:
        """Two-dimensional control-volume area in m²."""

        return self.dx * self.dy

    @property
    def x_face_area(self) -> float:
        """Length of a face normal to the x direction in m."""

        return self.dy

    @property
    def y_face_area(self) -> float:
        """Length of a face normal to the y direction in m."""

        return self.dx

    @property
    def x_centers(self) -> NDArray[np.float64]:
        """One-dimensional x coordinates of cell centers in m."""

        return self.x_min + (np.arange(self.nx, dtype=np.float64) + 0.5) * self.dx

    @property
    def y_centers(self) -> NDArray[np.float64]:
        """One-dimensional y coordinates of cell centers in m."""

        return self.y_min + (np.arange(self.ny, dtype=np.float64) + 0.5) * self.dy

    @property
    def x_faces(self) -> NDArray[np.float64]:
        """One-dimensional x coordinates of cell faces in m."""

        return np.linspace(self.x_min, self.x_max, self.nx + 1, dtype=np.float64)

    @property
    def y_faces(self) -> NDArray[np.float64]:
        """One-dimensional y coordinates of cell faces in m."""

        return np.linspace(self.y_min, self.y_max, self.ny + 1, dtype=np.float64)

    @property
    def cell_centers(
        self,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Cell-center coordinate arrays ``(x, y)`` with shape ``(ny, nx)``."""

        x, y = np.meshgrid(self.x_centers, self.y_centers, indexing="xy")
        return np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)

    @property
    def x_face_centers(
        self,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Coordinates of x-normal faces with shape ``(ny, nx + 1)``."""

        x, y = np.meshgrid(self.x_faces, self.y_centers, indexing="xy")
        return np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)

    @property
    def y_face_centers(
        self,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Coordinates of y-normal faces with shape ``(ny + 1, nx)``."""

        x, y = np.meshgrid(self.x_centers, self.y_faces, indexing="xy")
        return np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
