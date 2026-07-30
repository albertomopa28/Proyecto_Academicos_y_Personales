"""Explicit ghost-cell boundary conditions."""

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError

BoundaryCondition = Literal["periodic", "transmissive"]


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
