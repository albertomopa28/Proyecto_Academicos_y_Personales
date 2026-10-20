"""Navas-Montilla mixture, radiation, wind, and topographic relations."""

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError
from wildfire_model.mesh.mesh_2d import UniformMesh2D

type MixtureApproximation = Literal["complete", "initial", "fuel_only"]


def _compatible_fields(
    named_values: tuple[tuple[str, ArrayLike], ...],
) -> tuple[NDArray[np.float64], ...]:
    arrays = tuple(np.asarray(value, dtype=np.float64) for _, value in named_values)
    shapes = {array.shape for array in arrays if array.ndim > 0}
    if len(shapes) > 1:
        details = ", ".join(
            f"{name}={array.shape}" for (name, _), array in zip(named_values, arrays, strict=True)
        )
        raise InvalidConfigurationError(f"Incompatible field shapes: {details}.")
    target_shape = next(iter(shapes), ())
    expanded = tuple(
        np.full(target_shape, float(array), dtype=np.float64)
        if array.ndim == 0 and target_shape
        else array.copy()
        for array in arrays
    )
    if any(not np.all(np.isfinite(array)) for array in expanded):
        raise InvalidConfigurationError("Physical property fields must be finite.")
    return expanded


def initial_fuel_load(
    maximum_fuel_load: ArrayLike,
    surface_coverage: ArrayLike,
) -> NDArray[np.float64]:
    """Return Navas-Montilla Sec. 2.1, Eq. (4), in kg/m²."""

    maximum, coverage = _compatible_fields(
        (
            ("maximum fuel load", maximum_fuel_load),
            ("surface coverage", surface_coverage),
        )
    )
    if np.any(maximum < 0.0) or np.any((coverage < 0.0) | (coverage > 1.0)):
        raise InvalidConfigurationError(
            "Maximum fuel load must be non-negative and coverage must lie in [0, 1]."
        )
    return maximum * coverage


def initial_fuel_volume_fraction(
    fuel_load: ArrayLike,
    intrinsic_fuel_density: ArrayLike,
    layer_height: float,
) -> NDArray[np.float64]:
    """Return Navas-Montilla Sec. 2.1, Eq. (3)."""

    load, density = _compatible_fields(
        (("fuel load", fuel_load), ("intrinsic fuel density", intrinsic_fuel_density))
    )
    if (
        np.any(load < 0.0)
        or np.any(density <= 0.0)
        or not np.isfinite(layer_height)
        or layer_height <= 0.0
    ):
        raise InvalidConfigurationError(
            "Fuel load must be non-negative; density and layer height must be positive."
        )
    ratio = load / (density * layer_height)
    if np.any(ratio > 1.0):
        raise InvalidConfigurationError("Initial fuel volume fraction must not exceed one.")
    return ratio


def mixture_density(
    intrinsic_fuel_density: ArrayLike,
    air_density: ArrayLike,
    initial_fuel_fraction: ArrayLike,
    biomass: ArrayLike,
) -> NDArray[np.float64]:
    """Return Navas-Montilla Sec. 2.1, Eq. (6), in kg/m³."""

    fuel_density, air, initial_fraction, remaining = _compatible_fields(
        (
            ("intrinsic fuel density", intrinsic_fuel_density),
            ("air density", air_density),
            ("initial fuel fraction", initial_fuel_fraction),
            ("biomass", biomass),
        )
    )
    if (
        np.any(fuel_density <= 0.0)
        or np.any(air <= 0.0)
        or np.any((initial_fraction < 0.0) | (initial_fraction > 1.0))
        or np.any((remaining < 0.0) | (remaining > 1.0))
    ):
        raise InvalidConfigurationError("Mixture inputs violate physical bounds.")
    return (fuel_density - air) * initial_fraction * remaining + air


def mixture_specific_heat(
    intrinsic_fuel_density: ArrayLike,
    fuel_specific_heat: ArrayLike,
    air_density: ArrayLike,
    air_specific_heat: ArrayLike,
    initial_fuel_fraction: ArrayLike,
    biomass: ArrayLike,
) -> NDArray[np.float64]:
    """Return Navas-Montilla Sec. 2.1, Eq. (7), in J/(kg K)."""

    (
        fuel_density,
        fuel_heat,
        air,
        air_heat,
        initial_fraction,
        remaining,
    ) = _compatible_fields(
        (
            ("intrinsic fuel density", intrinsic_fuel_density),
            ("fuel specific heat", fuel_specific_heat),
            ("air density", air_density),
            ("air specific heat", air_specific_heat),
            ("initial fuel fraction", initial_fuel_fraction),
            ("biomass", biomass),
        )
    )
    density = mixture_density(fuel_density, air, initial_fraction, remaining)
    if np.any(fuel_heat <= 0.0) or np.any(air_heat <= 0.0):
        raise InvalidConfigurationError("Specific heats must be positive.")
    numerator = (
        fuel_heat * fuel_density - air_heat * air
    ) * initial_fraction * remaining + air_heat * air
    return numerator / density


def initial_mixture_properties(
    intrinsic_fuel_density: ArrayLike,
    fuel_specific_heat: ArrayLike,
    air_density: ArrayLike,
    air_specific_heat: ArrayLike,
    initial_fuel_fraction: ArrayLike,
    approximation: MixtureApproximation = "complete",
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return the Eq. (8) mixture or the explicit Eq. (9) approximation."""

    if approximation in ("complete", "initial"):
        density = mixture_density(
            intrinsic_fuel_density,
            air_density,
            initial_fuel_fraction,
            1.0,
        )
        specific_heat = mixture_specific_heat(
            intrinsic_fuel_density,
            fuel_specific_heat,
            air_density,
            air_specific_heat,
            initial_fuel_fraction,
            1.0,
        )
        return density, specific_heat
    if approximation == "fuel_only":
        fuel_density, fuel_heat, initial_fraction = _compatible_fields(
            (
                ("intrinsic fuel density", intrinsic_fuel_density),
                ("fuel specific heat", fuel_specific_heat),
                ("initial fuel fraction", initial_fuel_fraction),
            )
        )
        if (
            np.any(fuel_density <= 0.0)
            or np.any(fuel_heat <= 0.0)
            or np.any((initial_fraction < 0.0) | (initial_fraction > 1.0))
        ):
            raise InvalidConfigurationError("Fuel-only properties are nonphysical.")
        return fuel_density * initial_fraction, fuel_heat
    raise InvalidConfigurationError(f"Unsupported mixture approximation {approximation!r}.")


def mixture_properties(
    intrinsic_fuel_density: ArrayLike,
    fuel_specific_heat: ArrayLike,
    air_density: ArrayLike,
    air_specific_heat: ArrayLike,
    initial_fuel_fraction: ArrayLike,
    biomass: ArrayLike,
    approximation: MixtureApproximation = "complete",
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Select Eqs. (6)-(7), the initial Eq. (8), or simplified Eq. (9)."""

    if approximation == "complete":
        density = mixture_density(
            intrinsic_fuel_density,
            air_density,
            initial_fuel_fraction,
            biomass,
        )
        specific_heat = mixture_specific_heat(
            intrinsic_fuel_density,
            fuel_specific_heat,
            air_density,
            air_specific_heat,
            initial_fuel_fraction,
            biomass,
        )
        return density, specific_heat
    return initial_mixture_properties(
        intrinsic_fuel_density,
        fuel_specific_heat,
        air_density,
        air_specific_heat,
        initial_fuel_fraction,
        approximation,
    )


def effective_thermal_conductivity(
    temperature: ArrayLike,
    conductive_thermal_conductivity: ArrayLike,
    emissivity: float,
    optical_path_length: ArrayLike,
    stefan_boltzmann_constant: float = 5.670e-8,
) -> NDArray[np.float64]:
    """Return Navas-Montilla Sec. 2.4, Eq. (28), in W/(m K)."""

    values, conductive, path = _compatible_fields(
        (
            ("temperature", temperature),
            ("conductive thermal conductivity", conductive_thermal_conductivity),
            ("optical path length", optical_path_length),
        )
    )
    if (
        np.any(values <= 0.0)
        or np.any(conductive < 0.0)
        or np.any(path < 0.0)
        or not np.isfinite(emissivity)
        or not 0.0 <= emissivity <= 1.0
        or not np.isfinite(stefan_boltzmann_constant)
        or stefan_boltzmann_constant <= 0.0
    ):
        raise InvalidConfigurationError("Effective-conductivity inputs are nonphysical.")
    return conductive + 4.0 * stefan_boltzmann_constant * emissivity * path * values**3


def wind_attenuation_factor(
    fuel_fraction: ArrayLike,
    intrinsic_fuel_density: ArrayLike,
    fuel_specific_heat: ArrayLike,
    air_density: ArrayLike,
    air_specific_heat: ArrayLike,
) -> NDArray[np.float64]:
    """Return Navas-Montilla Sec. 2.3, Eq. (39)."""

    fraction, fuel_density, fuel_heat, air, air_heat = _compatible_fields(
        (
            ("fuel fraction", fuel_fraction),
            ("intrinsic fuel density", intrinsic_fuel_density),
            ("fuel specific heat", fuel_specific_heat),
            ("air density", air_density),
            ("air specific heat", air_specific_heat),
        )
    )
    if (
        np.any((fraction < 0.0) | (fraction > 1.0))
        or np.any(fuel_density <= 0.0)
        or np.any(fuel_heat <= 0.0)
        or np.any(air <= 0.0)
        or np.any(air_heat <= 0.0)
    ):
        raise InvalidConfigurationError("Wind-attenuation inputs are nonphysical.")
    air_capacity = (1.0 - fraction) * air * air_heat
    total_capacity = fraction * fuel_density * fuel_heat + air_capacity
    return air_capacity / total_capacity


def elevation_gradient(
    elevation: ArrayLike,
    mesh: UniformMesh2D,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return the dimensionless Cartesian gradient of elevation in m/m."""

    values = np.asarray(elevation, dtype=np.float64)
    if values.shape != mesh.shape or not np.all(np.isfinite(values)):
        raise InvalidConfigurationError(f"Elevation must be finite with shape {mesh.shape}.")
    gradients = []
    for axis, spacing in ((0, mesh.dy), (1, mesh.dx)):
        if values.shape[axis] == 1:
            gradients.append(np.zeros_like(values))
        else:
            gradients.append(
                np.gradient(
                    values, spacing, axis=axis, edge_order=2 if values.shape[axis] >= 3 else 1
                )
            )
    gradient_y, gradient_x = gradients
    return (
        np.asarray(gradient_x, dtype=np.float64),
        np.asarray(gradient_y, dtype=np.float64),
    )


def effective_velocity(
    wind_x: ArrayLike,
    wind_y: ArrayLike,
    elevation: ArrayLike,
    mesh: UniformMesh2D,
    *,
    wind_correction: float,
    slope_correction: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return Eq. (38), interpreting slope correction in m/s."""

    wind_x_array, wind_y_array = _compatible_fields((("x wind", wind_x), ("y wind", wind_y)))
    if wind_x_array.shape == ():
        wind_x_array = np.full(mesh.shape, float(wind_x_array))
        wind_y_array = np.full(mesh.shape, float(wind_y_array))
    if wind_x_array.shape != mesh.shape:
        raise InvalidConfigurationError(f"Wind fields must be scalar or have shape {mesh.shape}.")
    if (
        not np.isfinite(wind_correction)
        or wind_correction < 0.0
        or not np.isfinite(slope_correction)
        or slope_correction < 0.0
    ):
        raise InvalidConfigurationError("Velocity corrections must be non-negative.")
    gradient_x, gradient_y = elevation_gradient(elevation, mesh)
    return (
        wind_correction * wind_x_array + slope_correction * gradient_x,
        wind_correction * wind_y_array + slope_correction * gradient_y,
    )
