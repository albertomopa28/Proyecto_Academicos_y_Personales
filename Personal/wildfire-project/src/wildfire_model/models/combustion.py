"""Combustion-rate functions for the base Reisch model."""

from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError, NonPhysicalStateError
from wildfire_model.models.parameters import ADRParameters

type CombustionModel = Literal["arrhenius", "constant"]


def arrhenius_rate(
    temperature: ArrayLike,
    parameters: ADRParameters,
) -> NDArray[np.float64]:
    r"""Return the switched Arrhenius rate :math:`\Psi(T)` in s\ :sup:`-1`.

    This is Reisch et al. (2024), Sec. 2.1, Eqs. (9)--(10):

    ``Psi(T) = s(T) A exp(-T_ac/T)``, with ``s=0`` below ``T_pc`` and
    ``s=1`` at and above ``T_pc``. No smoothing or regularization is
    applied to the published discontinuous switch.
    """

    temperature_array = np.asarray(temperature, dtype=np.float64)
    if not np.all(np.isfinite(temperature_array)) or np.any(temperature_array <= 0.0):
        raise NonPhysicalStateError(
            "Arrhenius rate requires finite absolute temperatures strictly above 0 K."
        )
    active = temperature_array >= parameters.pyrolysis_temperature
    rate = np.zeros_like(temperature_array)
    rate[active] = parameters.pre_exponential_factor * np.exp(
        -parameters.activation_temperature / temperature_array[active]
    )
    return rate


def constant_combustion_rate(
    temperature: ArrayLike,
    pre_exponential_factor: float,
    pyrolysis_temperature: float,
) -> NDArray[np.float64]:
    """Return Navas-Montilla Sec. 2.5, Eqs. (41)-(42), in s⁻¹."""

    temperature_array = np.asarray(temperature, dtype=np.float64)
    if not np.all(np.isfinite(temperature_array)) or np.any(temperature_array <= 0.0):
        raise NonPhysicalStateError(
            "Combustion rate requires finite absolute temperatures above 0 K."
        )
    if not np.isfinite(pre_exponential_factor) or pre_exponential_factor < 0.0:
        raise InvalidConfigurationError(
            "Constant combustion pre-exponential factor must be finite and non-negative."
        )
    if not np.isfinite(pyrolysis_temperature) or pyrolysis_temperature <= 0.0:
        raise InvalidConfigurationError("Pyrolysis temperature must be finite and positive.")
    return np.where(
        temperature_array >= pyrolysis_temperature,
        pre_exponential_factor,
        0.0,
    )
