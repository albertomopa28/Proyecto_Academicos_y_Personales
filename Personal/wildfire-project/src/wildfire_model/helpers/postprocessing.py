"""Front-position, propagation-speed, error, and convergence diagnostics."""

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError
from wildfire_model.mesh.mesh_2d import UniformMesh2D

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


def _checked_2d_field(
    values: ArrayLike,
    mesh: UniformMesh2D,
    name: str,
) -> NDArray[np.float64]:
    array = np.asarray(values, dtype=np.float64)
    if array.shape != mesh.shape or not np.all(np.isfinite(array)):
        raise InvalidConfigurationError(f"{name} must be finite with shape {mesh.shape}.")
    return array


def _interface_mask(region: NDArray[np.bool_]) -> NDArray[np.bool_]:
    padded = np.pad(region, 1, mode="edge")
    interior = padded[1:-1, 1:-1]
    surrounded = padded[1:-1, :-2] & padded[1:-1, 2:] & padded[:-2, 1:-1] & padded[2:, 1:-1]
    return interior & ~surrounded


def isotherm_front(
    temperature: ArrayLike,
    mesh: UniformMesh2D,
    threshold: float,
) -> NDArray[np.bool_]:
    """Return hot-side cells adjacent to a declared isotherm."""

    values = _checked_2d_field(temperature, mesh, "temperature")
    if not np.isfinite(threshold):
        raise InvalidConfigurationError("Isotherm threshold must be finite.")
    return _interface_mask(values >= threshold)


def consumed_biomass_front(
    biomass: ArrayLike,
    mesh: UniformMesh2D,
    threshold: float = 0.99,
) -> NDArray[np.bool_]:
    """Return consumed-side cells adjacent to the biomass threshold."""

    values = _checked_2d_field(biomass, mesh, "biomass")
    if not np.isfinite(threshold) or not 0.0 <= threshold <= 1.0:
        raise InvalidConfigurationError("Biomass threshold must lie in [0, 1].")
    return _interface_mask(values < threshold)


def burned_area(
    biomass: ArrayLike,
    mesh: UniformMesh2D,
    threshold: float = 0.99,
) -> float:
    """Return area in m² with biomass below the declared threshold."""

    values = _checked_2d_field(biomass, mesh, "biomass")
    if not np.isfinite(threshold) or not 0.0 <= threshold <= 1.0:
        raise InvalidConfigurationError("Biomass threshold must lie in [0, 1].")
    return float(np.count_nonzero(values < threshold) * mesh.cell_area)


def front_extents(
    front: ArrayLike,
    mesh: UniformMesh2D,
) -> tuple[float, float, float, float]:
    """Return left, right, lower, and upper front coordinates in m."""

    mask = np.asarray(front, dtype=np.bool_)
    if mask.shape != mesh.shape or not np.any(mask):
        raise InvalidConfigurationError(
            f"Front mask must have shape {mesh.shape} and contain a front."
        )
    x, y = mesh.cell_centers
    return (
        float(np.min(x[mask])),
        float(np.max(x[mask])),
        float(np.min(y[mask])),
        float(np.max(y[mask])),
    )


def burned_centroid(
    biomass: ArrayLike,
    mesh: UniformMesh2D,
) -> tuple[float, float]:
    """Return the consumption-weighted centroid in m."""

    values = _checked_2d_field(biomass, mesh, "biomass")
    if np.any((values < 0.0) | (values > 1.0)):
        raise InvalidConfigurationError("Biomass must lie in [0, 1].")
    weights = 1.0 - values
    total = float(np.sum(weights))
    if total <= 0.0:
        raise InvalidConfigurationError("Burned centroid requires consumed biomass.")
    x, y = mesh.cell_centers
    return float(np.sum(weights * x) / total), float(np.sum(weights * y) / total)


def front_anisotropy(
    front: ArrayLike,
    mesh: UniformMesh2D,
) -> float:
    """Return x span divided by y span for a non-degenerate front."""

    left, right, lower, upper = front_extents(front, mesh)
    y_span = upper - lower
    if y_span <= 0.0:
        raise InvalidConfigurationError("Front anisotropy requires a positive y span.")
    return (right - left) / y_span
