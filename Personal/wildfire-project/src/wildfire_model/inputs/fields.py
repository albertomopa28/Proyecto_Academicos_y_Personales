"""Validation and evaluation of scalar and vector fields on a 2D mesh."""

from collections.abc import Callable
from inspect import signature
from math import isfinite
from typing import cast

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError
from wildfire_model.mesh.mesh_2d import UniformMesh2D

type SpatialField = Callable[[NDArray[np.float64], NDArray[np.float64]], ArrayLike]
type TimeDependentField = Callable[[NDArray[np.float64], NDArray[np.float64], float], ArrayLike]
type FieldInput = float | ArrayLike | SpatialField | TimeDependentField


def evaluate_field(
    field: FieldInput,
    mesh: UniformMesh2D,
    *,
    name: str,
    time: float | None = None,
    minimum: float | None = None,
    maximum: float | None = None,
) -> NDArray[np.float64]:
    """Evaluate a scalar, array, or callable as a finite ``(ny, nx)`` field.

    Scalars are expanded explicitly. Arrays must already have the canonical
    mesh shape; incompatible arrays are never broadcast.
    """

    if time is not None and not isfinite(time):
        raise InvalidConfigurationError(f"{name} evaluation time must be finite.")
    if callable(field):
        x, y = mesh.cell_centers
        parameter_count = len(signature(field).parameters)
        if parameter_count == 2:
            values = cast(SpatialField, field)(x, y)
        elif parameter_count == 3 and time is not None:
            values = cast(TimeDependentField, field)(x, y, time)
        elif parameter_count == 3:
            raise InvalidConfigurationError(
                f"{name} is time-dependent but no evaluation time was supplied."
            )
        else:
            raise InvalidConfigurationError(f"{name} callable must accept (x, y) or (x, y, time).")
    else:
        values = field
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 0:
        array = np.full(mesh.shape, float(array), dtype=np.float64)
    elif array.shape != mesh.shape:
        raise InvalidConfigurationError(
            f"{name} must be scalar or have shape {mesh.shape}; received {array.shape}."
        )
    else:
        array = array.copy()
    if not np.all(np.isfinite(array)):
        raise InvalidConfigurationError(f"{name} contains NaN or infinity.")
    if minimum is not None and np.any(array < minimum):
        raise InvalidConfigurationError(f"{name} must be >= {minimum:.16g}.")
    if maximum is not None and np.any(array > maximum):
        raise InvalidConfigurationError(f"{name} must be <= {maximum:.16g}.")
    return array


def evaluate_vector_field(
    x_component: FieldInput,
    y_component: FieldInput,
    mesh: UniformMesh2D,
    *,
    name: str,
    time: float | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Evaluate both components of a velocity field in m/s."""

    return (
        evaluate_field(x_component, mesh, name=f"{name} x component", time=time),
        evaluate_field(y_component, mesh, name=f"{name} y component", time=time),
    )
