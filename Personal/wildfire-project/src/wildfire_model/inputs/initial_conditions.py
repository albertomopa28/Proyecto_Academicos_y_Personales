"""Reusable two-dimensional temperature and biomass initial conditions."""

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError
from wildfire_model.inputs.fields import FieldInput, evaluate_field
from wildfire_model.mesh.mesh_2d import UniformMesh2D
from wildfire_model.models.state import ADRState


def circular_ignition(
    mesh: UniformMesh2D,
    center: tuple[float, float],
    radius: float,
    *,
    ignition_temperature: float,
    ambient_temperature: float,
    biomass: FieldInput = 1.0,
) -> ADRState:
    """Initialize a circular hot region, as in Navas-Montilla Eqs. (80)-(82)."""

    if not np.isfinite(radius) or radius <= 0.0:
        raise InvalidConfigurationError("Ignition radius must be finite and positive.")
    if (
        not np.isfinite(ignition_temperature)
        or not np.isfinite(ambient_temperature)
        or ignition_temperature <= 0.0
        or ambient_temperature <= 0.0
    ):
        raise InvalidConfigurationError(
            "Ignition and ambient temperatures must be finite and positive."
        )
    x, y = mesh.cell_centers
    distance = np.hypot(x - center[0], y - center[1])
    temperature = np.where(distance < radius, ignition_temperature, ambient_temperature)
    biomass_field = evaluate_field(
        biomass,
        mesh,
        name="initial biomass",
        minimum=0.0,
        maximum=1.0,
    )
    return ADRState.from_arrays(temperature, biomass_field)


def gaussian_ignition(
    mesh: UniformMesh2D,
    center: tuple[float, float],
    width: float,
    *,
    peak_temperature: float,
    ambient_temperature: float,
    biomass: FieldInput = 1.0,
) -> ADRState:
    """Initialize a Gaussian temperature perturbation with width in m."""

    if not np.isfinite(width) or width <= 0.0:
        raise InvalidConfigurationError("Gaussian ignition width must be positive.")
    if peak_temperature <= ambient_temperature or ambient_temperature <= 0.0:
        raise InvalidConfigurationError(
            "Peak temperature must exceed the positive ambient temperature."
        )
    x, y = mesh.cell_centers
    radius_squared = (x - center[0]) ** 2 + (y - center[1]) ** 2
    temperature = ambient_temperature + (peak_temperature - ambient_temperature) * np.exp(
        -0.5 * radius_squared / width**2
    )
    biomass_field = evaluate_field(
        biomass,
        mesh,
        name="initial biomass",
        minimum=0.0,
        maximum=1.0,
    )
    return ADRState.from_arrays(temperature, biomass_field)


def line_ignition(
    mesh: UniformMesh2D,
    start: tuple[float, float],
    end: tuple[float, float],
    half_width: float,
    *,
    ignition_temperature: float,
    ambient_temperature: float,
    biomass: FieldInput = 1.0,
) -> ADRState:
    """Initialize a finite ignition line with a declared half-width in m."""

    if not np.isfinite(half_width) or half_width <= 0.0:
        raise InvalidConfigurationError("Ignition-line half-width must be positive.")
    direction = np.asarray(end, dtype=np.float64) - np.asarray(start, dtype=np.float64)
    length_squared = float(direction @ direction)
    if length_squared == 0.0 or not np.isfinite(length_squared):
        raise InvalidConfigurationError("Ignition-line endpoints must be finite and distinct.")
    x, y = mesh.cell_centers
    relative_x = x - start[0]
    relative_y = y - start[1]
    projection = np.clip(
        (relative_x * direction[0] + relative_y * direction[1]) / length_squared,
        0.0,
        1.0,
    )
    closest_x = start[0] + projection * direction[0]
    closest_y = start[1] + projection * direction[1]
    hot = np.hypot(x - closest_x, y - closest_y) <= half_width
    temperature = np.where(hot, ignition_temperature, ambient_temperature)
    biomass_field = evaluate_field(
        biomass,
        mesh,
        name="initial biomass",
        minimum=0.0,
        maximum=1.0,
    )
    return ADRState.from_arrays(temperature, biomass_field)


def apply_firebreak(
    biomass: ArrayLike,
    mesh: UniformMesh2D,
    mask: FieldInput,
    *,
    remaining_fraction: float = 0.0,
) -> NDArray[np.float64]:
    """Set masked fuel to a zero or reduced remaining fraction."""

    biomass_field = evaluate_field(
        biomass,
        mesh,
        name="biomass",
        minimum=0.0,
        maximum=1.0,
    )
    mask_field = evaluate_field(mask, mesh, name="firebreak mask", minimum=0.0)
    if not np.isfinite(remaining_fraction) or not 0.0 <= remaining_fraction <= 1.0:
        raise InvalidConfigurationError("remaining_fraction must lie in [0, 1].")
    return np.where(mask_field.astype(np.bool_), remaining_fraction, biomass_field)
