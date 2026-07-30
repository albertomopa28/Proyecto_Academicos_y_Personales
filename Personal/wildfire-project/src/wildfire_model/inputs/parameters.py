"""Predefined physical parameter sets."""

from wildfire_model.models.parameters import ADRParameters


def reisch_reaction_parameters() -> ADRParameters:
    """Return Reisch et al. Table 1 defaults for Appendix A.1."""

    return ADRParameters()


def reisch_case_2_parameters() -> ADRParameters:
    """Return Reisch Appendix A.2 advection-diffusion-cooling parameters."""

    return ADRParameters(
        bulk_density=1.0,
        specific_heat=1_000.0,
        pre_exponential_factor=0.0,
        heating_value=4_000_000.0,
        cooling_coefficient=10.0,
        ambient_temperature=300.0,
        pyrolysis_temperature=400.0,
        activation_temperature=400.0,
        thermal_conductivity=10_000.0,
        velocity=5.0,
    )
