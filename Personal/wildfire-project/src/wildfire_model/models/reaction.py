"""Spatially homogeneous temperature--biomass reaction system.

The local system is Reisch et al. (2024), Sec. 3.1, Eq. (16), obtained
from the reaction and ambient-cooling terms of Eq. (14):

``rho_0 C dT/dt = -h(T-T_inf) + Psi(T) rho_0 H Y``
``dY/dt = -Psi(T)Y``.

"""

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.models.combustion import arrhenius_rate
from wildfire_model.models.parameters import ADRParameters
from wildfire_model.models.state import validate_temperature_biomass


def reaction_rhs(
    temperature: ArrayLike,
    biomass: ArrayLike,
    parameters: ADRParameters,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Evaluate the dimensional right-hand side ``(dT/dt, dY/dt)``."""

    temperature_array, biomass_array = validate_temperature_biomass(
        temperature, biomass, context="reaction right-hand side"
    )
    rate = arrhenius_rate(temperature_array, parameters)
    biomass_rate = -rate * biomass_array
    temperature_rate = (
        -parameters.cooling_rate * (temperature_array - parameters.ambient_temperature)
        - parameters.heat_release_temperature * biomass_rate
    )
    return temperature_rate, biomass_rate
