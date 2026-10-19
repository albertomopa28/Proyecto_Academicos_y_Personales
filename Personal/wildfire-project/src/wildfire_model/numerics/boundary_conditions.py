"""Explicit ghost-cell boundary conditions."""

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError

BoundaryCondition = Literal["periodic", "transmissive"]
BoundaryCondition2D = Literal[
    "periodic",
    "transmissive",
    "zero_diffusive_flux",
    "zero_total_flux",
]


def extend_with_ghost_cells(
    values: ArrayLike,
    ghost_width: int,
    boundary: BoundaryCondition,
) -> NDArray[np.float64]:
    """Return cell averages with periodic or zero-gradient ghost cells.

    ``transmissive`` copies the nearest interior cell into every ghost
    cell. It is therefore a discrete zero-normal-gradient condition, not
    the zero-total-heat-flux condition in Reisch et al. Eq. (15).
    """

    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size == 0:
        raise InvalidConfigurationError("Ghost-cell input must be a non-empty 1D array.")
    if not np.all(np.isfinite(array)):
        raise InvalidConfigurationError("Ghost-cell input contains NaN or infinity.")
    if isinstance(ghost_width, bool) or ghost_width < 0:
        raise InvalidConfigurationError("ghost_width must be a non-negative integer.")
    if ghost_width == 0:
        return array.copy()
    if boundary == "periodic":
        if ghost_width > array.size:
            repeats = int(np.ceil(ghost_width / array.size))
            tiled = np.tile(array, 2 * repeats + 1)
            middle = repeats * array.size
            return tiled[middle - ghost_width : middle + array.size + ghost_width].copy()
        return np.concatenate((array[-ghost_width:], array, array[:ghost_width]))
    if boundary == "transmissive":
        return np.pad(array, ghost_width, mode="edge")
    raise InvalidConfigurationError(
        f"Unsupported boundary condition {boundary!r}; expected 'periodic' or 'transmissive'."
    )


def extend_with_ghost_cells_2d(
    values: ArrayLike,
    ghost_width: int,
    boundary: BoundaryCondition2D,
) -> NDArray[np.float64]:
    """Return a 2D field with ghost cells on both coordinate axes.

    ``transmissive`` and ``zero_diffusive_flux`` both use zero-gradient
    ghost values, but retain distinct boundary names because transmissive
    advection is not a zero-total-flux condition. ``zero_total_flux`` is
    enforced on the combined finite-volume flux after reconstruction.
    """

    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 2 or array.size == 0:
        raise InvalidConfigurationError("Ghost-cell input must be a non-empty 2D array.")
    if not np.all(np.isfinite(array)):
        raise InvalidConfigurationError("Ghost-cell input contains NaN or infinity.")
    if isinstance(ghost_width, bool) or ghost_width < 0:
        raise InvalidConfigurationError("ghost_width must be a non-negative integer.")
    if ghost_width == 0:
        return array.copy()
    if boundary == "periodic":
        return np.pad(array, ghost_width, mode="wrap")
    if boundary in ("transmissive", "zero_diffusive_flux", "zero_total_flux"):
        return np.pad(array, ghost_width, mode="edge")
    raise InvalidConfigurationError(f"Unsupported 2D boundary condition {boundary!r}.")


def reconstruction_boundary(boundary: BoundaryCondition2D) -> BoundaryCondition:
    """Map physical-flux boundary names to their ghost-state treatment."""

    return "periodic" if boundary == "periodic" else "transmissive"
