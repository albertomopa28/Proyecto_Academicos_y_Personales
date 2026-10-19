"""Navas-Montilla apparent heat-capacity moisture relations."""

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError

type MoistureModel = Literal["NM", "S2M", "S3M"]


def _fields(*values: ArrayLike) -> tuple[NDArray[np.float64], ...]:
    arrays = tuple(np.asarray(value, dtype=np.float64) for value in values)
    shapes = {array.shape for array in arrays if array.ndim > 0}
    if len(shapes) > 1:
        raise InvalidConfigurationError("Moisture fields have incompatible shapes.")
    shape = next(iter(shapes), ())
    result = tuple(
        np.full(shape, float(array), dtype=np.float64)
        if array.ndim == 0 and shape
        else array.copy()
        for array in arrays
    )
    if any(not np.all(np.isfinite(array)) for array in result):
        raise InvalidConfigurationError("Moisture inputs must be finite.")
    return result


def s2m_enthalpy_increments(
    moisture: ArrayLike,
    *,
    dry_specific_heat: float,
    water_specific_heat: float,
    latent_heat: float,
    ambient_temperature: float,
    evaporation_temperature: float,
    pyrolysis_temperature: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Return sensible, latent, and post-evaporation increments, Eqs. (45)-(47)."""

    moisture_array = np.asarray(moisture, dtype=np.float64)
    if (
        not np.all(np.isfinite(moisture_array))
        or np.any(moisture_array < 0.0)
        or min(dry_specific_heat, water_specific_heat, latent_heat) <= 0.0
        or not ambient_temperature < evaporation_temperature < pyrolysis_temperature
    ):
        raise InvalidConfigurationError("S2M inputs are nonphysical.")
    first = (dry_specific_heat + moisture_array * water_specific_heat) * (
        evaporation_temperature - ambient_temperature
    )
    latent = moisture_array * latent_heat
    second = np.full(
        moisture_array.shape,
        dry_specific_heat * (pyrolysis_temperature - evaporation_temperature),
        dtype=np.float64,
    )
    return first, latent, second


def s2m_effective_specific_heat(
    moisture: ArrayLike,
    *,
    dry_specific_heat: float,
    water_specific_heat: float,
    latent_heat: float,
    ambient_temperature: float,
    evaporation_temperature: float,
    pyrolysis_temperature: float,
) -> NDArray[np.float64]:
    """Return Navas-Montilla Sec. 2.6.2, Eq. (51), in J/(kg K)."""

    first, latent, second = s2m_enthalpy_increments(
        moisture,
        dry_specific_heat=dry_specific_heat,
        water_specific_heat=water_specific_heat,
        latent_heat=latent_heat,
        ambient_temperature=ambient_temperature,
        evaporation_temperature=evaporation_temperature,
        pyrolysis_temperature=pyrolysis_temperature,
    )
    return (first + latent + second) / (pyrolysis_temperature - ambient_temperature)


def s2m_apparent_specific_heat(
    temperature: ArrayLike,
    biomass: ArrayLike,
    moisture: ArrayLike,
    *,
    dry_specific_heat: float,
    water_specific_heat: float,
    latent_heat: float,
    ambient_temperature: float,
    evaporation_temperature: float,
    pyrolysis_temperature: float,
) -> NDArray[np.float64]:
    """Return the one-way S2M apparent capacity from Eq. (52)."""

    values, remaining, moisture_array = _fields(temperature, biomass, moisture)
    if (
        np.any(values <= 0.0)
        or np.any((remaining < 0.0) | (remaining > 1.0))
        or np.any(moisture_array < 0.0)
    ):
        raise InvalidConfigurationError("S2M state is nonphysical.")
    effective = s2m_effective_specific_heat(
        moisture_array,
        dry_specific_heat=dry_specific_heat,
        water_specific_heat=water_specific_heat,
        latent_heat=latent_heat,
        ambient_temperature=ambient_temperature,
        evaporation_temperature=evaporation_temperature,
        pyrolysis_temperature=pyrolysis_temperature,
    )
    return np.where(
        (values < pyrolysis_temperature) & (remaining == 1.0),
        effective,
        dry_specific_heat,
    )


def s3m_stage_specific_heats(
    moisture: ArrayLike,
    *,
    dry_specific_heat: float,
    water_specific_heat: float,
    latent_heat: float,
    ambient_temperature: float,
    evaporation_temperature: float,
    evaporation_interval: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return sensible and latent S3M capacities from Eqs. (53)-(57)."""

    moisture_array = np.asarray(moisture, dtype=np.float64)
    lower = evaporation_temperature - 0.5 * evaporation_interval
    upper = evaporation_temperature + 0.5 * evaporation_interval
    if (
        not np.all(np.isfinite(moisture_array))
        or np.any(moisture_array < 0.0)
        or not np.isfinite(evaporation_interval)
        or evaporation_interval <= 0.0
        or lower <= ambient_temperature
    ):
        raise InvalidConfigurationError("S3M inputs are nonphysical.")
    sensible = dry_specific_heat + moisture_array * water_specific_heat
    latent = (
        dry_specific_heat * (upper - lower)
        + moisture_array * (water_specific_heat * (evaporation_temperature - lower) + latent_heat)
    ) / (upper - lower)
    return sensible, latent


def s3m_apparent_specific_heat(
    temperature: ArrayLike,
    biomass: ArrayLike,
    moisture: ArrayLike,
    *,
    dry_specific_heat: float,
    water_specific_heat: float,
    latent_heat: float,
    ambient_temperature: float,
    evaporation_temperature: float,
    evaporation_interval: float,
) -> NDArray[np.float64]:
    """Return the one-way S3M apparent capacity from Eq. (58)."""

    values, remaining, moisture_array = _fields(temperature, biomass, moisture)
    if (
        np.any(values <= 0.0)
        or np.any((remaining < 0.0) | (remaining > 1.0))
        or np.any(moisture_array < 0.0)
    ):
        raise InvalidConfigurationError("S3M state is nonphysical.")
    sensible, latent = s3m_stage_specific_heats(
        moisture_array,
        dry_specific_heat=dry_specific_heat,
        water_specific_heat=water_specific_heat,
        latent_heat=latent_heat,
        ambient_temperature=ambient_temperature,
        evaporation_temperature=evaporation_temperature,
        evaporation_interval=evaporation_interval,
    )
    lower = evaporation_temperature - 0.5 * evaporation_interval
    upper = evaporation_temperature + 0.5 * evaporation_interval
    unburned = remaining == 1.0
    return np.where(
        unburned & (values <= lower),
        sensible,
        np.where(
            unburned & (values <= upper),
            latent,
            dry_specific_heat,
        ),
    )


def s3m_enthalpy(
    temperature: ArrayLike,
    moisture: ArrayLike,
    *,
    dry_specific_heat: float,
    water_specific_heat: float,
    latent_heat: float,
    ambient_temperature: float,
    evaporation_temperature: float,
    evaporation_interval: float,
) -> NDArray[np.float64]:
    """Integrate the S3M apparent capacity from ambient temperature."""

    values, moisture_array = _fields(temperature, moisture)
    sensible, latent = s3m_stage_specific_heats(
        moisture_array,
        dry_specific_heat=dry_specific_heat,
        water_specific_heat=water_specific_heat,
        latent_heat=latent_heat,
        ambient_temperature=ambient_temperature,
        evaporation_temperature=evaporation_temperature,
        evaporation_interval=evaporation_interval,
    )
    lower = evaporation_temperature - 0.5 * evaporation_interval
    upper = evaporation_temperature + 0.5 * evaporation_interval
    first = sensible * (lower - ambient_temperature)
    second = first + latent * (upper - lower)
    return np.where(
        values <= lower,
        sensible * (values - ambient_temperature),
        np.where(
            values <= upper,
            first + latent * (values - lower),
            second + dry_specific_heat * (values - upper),
        ),
    )


def dead_fuel_equilibrium_moisture(
    relative_humidity: ArrayLike,
    temperature_kelvin: ArrayLike,
) -> NDArray[np.float64]:
    """Return Navas-Montilla Sec. 2.6.4, Eqs. (61)-(62)."""

    humidity, temperature = _fields(relative_humidity, temperature_kelvin)
    if np.any((humidity < 0.0) | (humidity > 1.0)) or np.any(temperature <= 0.0):
        raise InvalidConfigurationError(
            "Relative humidity must lie in [0, 1] and temperature must be positive."
        )
    celsius = temperature - 273.15
    molecular_weight = 349.0 + 1.29 * celsius + 0.0135 * celsius**2
    k_value = 0.805 + 0.000736 * celsius - 0.00000273 * celsius**2
    k_one = 6.27 - 0.00938 * celsius - 0.000303 * celsius**2
    k_two = 1.91 + 0.04007 * celsius - 0.000293 * celsius**2
    product = k_value * humidity
    denominator_one = 1.0 - product
    denominator_two = 1.0 + k_one * product + k_one * k_two * product**2
    if np.any(denominator_one <= 0.0) or np.any(denominator_two <= 0.0):
        raise InvalidConfigurationError(
            "Dead-fuel equilibrium relation has a non-positive denominator."
        )
    return (
        18.0
        / molecular_weight
        * (
            product / denominator_one
            + (k_one * product + 2.0 * k_one * k_two * product**2) / denominator_two
        )
    )


def live_dead_moisture(
    live_fraction: ArrayLike,
    live_moisture: ArrayLike,
    dead_moisture: ArrayLike,
) -> NDArray[np.float64]:
    """Return the dry-basis live/dead mixture from Eqs. (59)-(60)."""

    fraction, live, dead = _fields(live_fraction, live_moisture, dead_moisture)
    if np.any((fraction < 0.0) | (fraction > 1.0)) or np.any(live < 0.0) or np.any(dead < 0.0):
        raise InvalidConfigurationError("Live/dead moisture inputs are nonphysical.")
    return fraction * live + (1.0 - fraction) * dead


def apparent_specific_heat(
    model: str,
    temperature: ArrayLike,
    biomass: ArrayLike,
    moisture: ArrayLike,
    *,
    dry_specific_heat: float,
    water_specific_heat: float,
    latent_heat: float,
    ambient_temperature: float,
    evaporation_temperature: float,
    pyrolysis_temperature: float,
    evaporation_interval: float | None = None,
) -> NDArray[np.float64]:
    """Evaluate a supported Navas-Montilla moisture model."""

    values, remaining, moisture_array = _fields(temperature, biomass, moisture)
    if model == "NM":
        return np.full(values.shape, dry_specific_heat, dtype=np.float64)
    if model == "S2M":
        return s2m_apparent_specific_heat(
            values,
            remaining,
            moisture_array,
            dry_specific_heat=dry_specific_heat,
            water_specific_heat=water_specific_heat,
            latent_heat=latent_heat,
            ambient_temperature=ambient_temperature,
            evaporation_temperature=evaporation_temperature,
            pyrolysis_temperature=pyrolysis_temperature,
        )
    if model == "S3M":
        if evaporation_interval is None:
            raise InvalidConfigurationError("S3M requires an explicit evaporation_interval in K.")
        return s3m_apparent_specific_heat(
            values,
            remaining,
            moisture_array,
            dry_specific_heat=dry_specific_heat,
            water_specific_heat=water_specific_heat,
            latent_heat=latent_heat,
            ambient_temperature=ambient_temperature,
            evaporation_temperature=evaporation_temperature,
            evaporation_interval=evaporation_interval,
        )
    if model == "C2M":
        raise InvalidConfigurationError(
            "C2M heat capacity is not exposed because Eq. (63) does not define "
            "the post-evaporation c_p' value or the SI units of Eq. (67)."
        )
    raise InvalidConfigurationError(f"Unsupported moisture model {model!r}.")
