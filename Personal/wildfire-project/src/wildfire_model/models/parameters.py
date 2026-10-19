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


@dataclass(frozen=True, slots=True)
class FuelProperties:
    """Navas-Montilla Tables 2 and 4 fuel properties in SI units."""

    intrinsic_density: float = 400.0
    dry_specific_heat: float = 1_000.0
    maximum_fuel_load: float = 4.0
    layer_height: float = 1.0
    pyrolysis_temperature: float = 600.0
    heating_value: float = 4_000_000.0

    def __post_init__(self) -> None:
        values = (
            self.intrinsic_density,
            self.dry_specific_heat,
            self.maximum_fuel_load,
            self.layer_height,
            self.pyrolysis_temperature,
            self.heating_value,
        )
        if not all(isfinite(value) and value > 0.0 for value in values):
            raise InvalidConfigurationError("Fuel properties must be finite and strictly positive.")


@dataclass(frozen=True, slots=True)
class EnvironmentalProperties:
    """Navas-Montilla ambient properties and fixed physical constants."""

    ambient_temperature: float = 300.0
    air_density: float = 1.0
    air_specific_heat: float = 1_000.0
    water_evaporation_temperature: float = 373.0
    water_latent_heat: float = 2_257_000.0
    water_specific_heat: float = 4_190.0
    stefan_boltzmann_constant: float = 5.670e-8
    emissivity: float = 0.9

    def __post_init__(self) -> None:
        values = (
            self.ambient_temperature,
            self.air_density,
            self.air_specific_heat,
            self.water_evaporation_temperature,
            self.water_latent_heat,
            self.water_specific_heat,
            self.stefan_boltzmann_constant,
            self.emissivity,
        )
        if not all(isfinite(value) and value > 0.0 for value in values):
            raise InvalidConfigurationError(
                "Environmental properties must be finite and strictly positive."
            )
        if self.emissivity > 1.0:
            raise InvalidConfigurationError("emissivity must not exceed one.")


@dataclass(frozen=True, slots=True)
class CalibrationParameters:
    """Navas-Montilla Table 3 calibration parameters in SI units."""

    conductive_thermal_conductivity: float = 0.1
    optical_path_length: float = 1.0
    surface_cooling_coefficient: float = 50.0
    wind_correction: float = 0.02
    slope_correction: float = 0.04

    def __post_init__(self) -> None:
        values = (
            self.conductive_thermal_conductivity,
            self.optical_path_length,
            self.surface_cooling_coefficient,
            self.wind_correction,
            self.slope_correction,
        )
        if not all(isfinite(value) and value >= 0.0 for value in values):
            raise InvalidConfigurationError(
                "Calibration parameters must be finite and non-negative."
            )
