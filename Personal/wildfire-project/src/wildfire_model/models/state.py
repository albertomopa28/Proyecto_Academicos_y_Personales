"""Physical-state validation shared by local and spatial models."""

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import NonPhysicalStateError


@dataclass(frozen=True, slots=True)
class ADRState:
    """Cell-average temperature [K] and remaining biomass fraction [-]."""

    temperature: NDArray[np.float64]
    biomass: NDArray[np.float64]

    @classmethod
    def from_arrays(cls, temperature: ArrayLike, biomass: ArrayLike) -> "ADRState":
        """Validate inputs and return an owned state."""

        temperature_array, biomass_array = validate_temperature_biomass(
            temperature, biomass, context="ADR state"
        )
        return cls(temperature_array.copy(), biomass_array.copy())


def validate_temperature_biomass(
    temperature: ArrayLike,
    biomass: ArrayLike,
    *,
    context: str,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Validate ``T > 0`` K and ``0 <= Y <= 1`` without altering the data."""

    temperature_array = np.asarray(temperature, dtype=np.float64)
    biomass_array = np.asarray(biomass, dtype=np.float64)
    if temperature_array.shape != biomass_array.shape:
        raise NonPhysicalStateError(
            f"{context}: temperature shape {temperature_array.shape} does not match "
            f"biomass shape {biomass_array.shape}."
        )
    if not np.all(np.isfinite(temperature_array)):
        raise NonPhysicalStateError(f"{context}: temperature contains NaN or infinity.")
    if not np.all(np.isfinite(biomass_array)):
        raise NonPhysicalStateError(f"{context}: biomass contains NaN or infinity.")
    if np.any(temperature_array <= 0.0):
        minimum = float(np.min(temperature_array))
        raise NonPhysicalStateError(
            f"{context}: absolute temperature must be > 0 K; minimum is {minimum:.16g} K."
        )
    if np.any((biomass_array < 0.0) | (biomass_array > 1.0)):
        minimum = float(np.min(biomass_array))
        maximum = float(np.max(biomass_array))
        raise NonPhysicalStateError(
            f"{context}: biomass must lie in [0, 1]; range is [{minimum:.16g}, {maximum:.16g}]."
        )
    return temperature_array, biomass_array
