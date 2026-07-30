"""Front-position, propagation-speed, error, and convergence diagnostics."""

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError

FrontDirection = Literal["left", "right"]


def front_position(
    coordinates: ArrayLike,
    temperature: ArrayLike,
    threshold: float,
    direction: FrontDirection = "right",
) -> float:
    """Interpolate the outermost crossing of a temperature threshold."""

    x = np.asarray(coordinates, dtype=np.float64)
    values = np.asarray(temperature, dtype=np.float64)
    if (
        x.ndim != 1
        or values.shape != x.shape
        or x.size < 2
        or np.any(~np.isfinite(x))
        or np.any(~np.isfinite(values))
        or np.any(np.diff(x) <= 0.0)
        or not np.isfinite(threshold)
    ):
        raise InvalidConfigurationError(
            "Front inputs must be finite 1D arrays on strictly increasing coordinates."
        )
    differences = values - threshold
    crossing_indices = np.flatnonzero(differences[:-1] * differences[1:] <= 0.0)
    crossing_indices = crossing_indices[
        differences[crossing_indices] != differences[crossing_indices + 1]
    ]
    if crossing_indices.size == 0:
        raise InvalidConfigurationError(
            f"Temperature profile has no crossing of threshold {threshold:.16g} K."
        )
    positions = x[crossing_indices] - differences[crossing_indices] * (
        x[crossing_indices + 1] - x[crossing_indices]
    ) / (differences[crossing_indices + 1] - differences[crossing_indices])
    if direction == "right":
        return float(np.max(positions))
    if direction == "left":
        return float(np.min(positions))
    raise InvalidConfigurationError("direction must be 'left' or 'right'.")


def propagation_speed(times: ArrayLike, positions: ArrayLike) -> float:
    """Return the least-squares front speed in m/s."""

    time_array = np.asarray(times, dtype=np.float64)
    position_array = np.asarray(positions, dtype=np.float64)
    if (
        time_array.ndim != 1
        or position_array.shape != time_array.shape
        or time_array.size < 2
        or np.any(~np.isfinite(time_array))
        or np.any(~np.isfinite(position_array))
        or np.any(np.diff(time_array) <= 0.0)
    ):
        raise InvalidConfigurationError(
            "Speed inputs must be finite equal-length arrays with increasing times."
        )
    design = np.column_stack((time_array, np.ones_like(time_array)))
    slope, _ = np.linalg.lstsq(design, position_array, rcond=None)[0]
    return float(slope)


def error_norms(
    numerical: ArrayLike,
    reference: ArrayLike,
) -> dict[str, float]:
    """Return discrete L1, L2, and Linfinity norms."""

    numerical_array = np.asarray(numerical, dtype=np.float64)
    reference_array = np.asarray(reference, dtype=np.float64)
    if numerical_array.shape != reference_array.shape or numerical_array.size == 0:
        raise InvalidConfigurationError(
            "Error-norm arrays must be non-empty and have identical shapes."
        )
    difference = numerical_array - reference_array
    return {
        "l1": float(np.mean(np.abs(difference))),
        "l2": float(np.sqrt(np.mean(difference**2))),
        "linf": float(np.max(np.abs(difference))),
    }


def observed_orders(
    resolutions: ArrayLike,
    errors: ArrayLike,
) -> NDArray[np.float64]:
    """Return pairwise ``log(e_i/e_{i+1}) / log(h_i/h_{i+1})``."""

    resolution_array = np.asarray(resolutions, dtype=np.float64)
    error_array = np.asarray(errors, dtype=np.float64)
    if (
        resolution_array.ndim != 1
        or error_array.ndim != 1
        or resolution_array.size != error_array.size
        or resolution_array.size < 2
        or np.any(~np.isfinite(resolution_array))
        or np.any(~np.isfinite(error_array))
        or np.any(resolution_array <= 0.0)
        or np.any(error_array <= 0.0)
        or np.any(np.diff(resolution_array) >= 0.0)
    ):
        raise InvalidConfigurationError(
            "resolutions and errors must be positive finite 1D arrays of equal length; "
            "resolutions must decrease."
        )
    return np.log(error_array[:-1] / error_array[1:]) / np.log(
        resolution_array[:-1] / resolution_array[1:]
    )
