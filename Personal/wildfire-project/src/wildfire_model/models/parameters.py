"""Dimensional parameters for the base Reisch ADR model.

Source
------
Reisch, Navas-Montilla & Özgen-Xian (2024), Sec. 2.1, Eqs. (9)--(14),
and Table 1. The paper tabulates energy in kJ and power in kW. This
module stores the exactly equivalent values in J and W so that all
internal quantities are SI.
"""

from dataclasses import dataclass
from math import isfinite

from wildfire_model.exceptions import InvalidConfigurationError


@dataclass(frozen=True, slots=True)
class ADRParameters:
    """Physical parameters for the one-dimensional base ADR system.

    Attributes use SI units. ``bulk_density`` is kg m^-3,
    ``specific_heat`` is J kg^-1 K^-1, ``pre_exponential_factor`` is
    s^-1, ``heating_value`` is J kg^-1, ``cooling_coefficient`` is
    W m^-3 K^-1, temperatures are K, ``thermal_conductivity`` is
    W m^-1 K^-1, and ``velocity`` is m s^-1.
    """

    bulk_density: float = 40.0
    specific_heat: float = 1_000.0
    pre_exponential_factor: float = 0.05
    heating_value: float = 4_000_000.0
    cooling_coefficient: float = 4_000.0
    ambient_temperature: float = 300.0
    pyrolysis_temperature: float = 400.0
    activation_temperature: float = 400.0
    thermal_conductivity: float = 2_000.0
    velocity: float = 0.0

    def __post_init__(self) -> None:
        values = {
            field: getattr(self, field)
            for field in (
                "bulk_density",
                "specific_heat",
                "pre_exponential_factor",
                "heating_value",
                "cooling_coefficient",
                "ambient_temperature",
                "pyrolysis_temperature",
                "activation_temperature",
                "thermal_conductivity",
                "velocity",
            )
        }
        nonfinite = [name for name, value in values.items() if not isfinite(value)]
        if nonfinite:
            raise InvalidConfigurationError(
                f"Physical parameters must be finite; invalid: {', '.join(nonfinite)}."
            )
        positive = {
            "bulk_density": self.bulk_density,
            "specific_heat": self.specific_heat,
            "heating_value": self.heating_value,
            "ambient_temperature": self.ambient_temperature,
            "pyrolysis_temperature": self.pyrolysis_temperature,
            "activation_temperature": self.activation_temperature,
        }
        invalid_positive = [name for name, value in positive.items() if value <= 0.0]
        if invalid_positive:
            raise InvalidConfigurationError(
                "Physical parameters must be strictly positive: "
                + ", ".join(invalid_positive)
                + "."
            )
        nonnegative = {
            "pre_exponential_factor": self.pre_exponential_factor,
            "cooling_coefficient": self.cooling_coefficient,
            "thermal_conductivity": self.thermal_conductivity,
        }
        invalid_nonnegative = [name for name, value in nonnegative.items() if value < 0.0]
        if invalid_nonnegative:
            raise InvalidConfigurationError(
                "Physical parameters cannot be negative: " + ", ".join(invalid_nonnegative) + "."
            )
        if self.ambient_temperature >= self.pyrolysis_temperature:
            raise InvalidConfigurationError(
                "ambient_temperature must be below pyrolysis_temperature for the "
                "Reisch ignition-switch model."
            )

    @property
    def volumetric_heat_capacity(self) -> float:
        """Return ``rho_0 C`` in J m^-3 K^-1."""

        return self.bulk_density * self.specific_heat

    @property
    def thermal_diffusivity(self) -> float:
        """Return ``k / (rho_0 C)`` in m^2 s^-1."""

        return self.thermal_conductivity / self.volumetric_heat_capacity

    @property
    def cooling_rate(self) -> float:
        """Return ``h / (rho_0 C)`` in s^-1."""

        return self.cooling_coefficient / self.volumetric_heat_capacity

    @property
    def heat_release_temperature(self) -> float:
        """Return ``H / C`` in K."""

        return self.heating_value / self.specific_heat
