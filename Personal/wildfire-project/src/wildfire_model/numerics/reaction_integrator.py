"""Implicit midpoint integration of the local reaction system.

The paper states that an implicit second-order Runge--Kutta method is
used (Sec. 2.2) but does not publish a Butcher tableau. This module
uses the one-stage Gauss method (implicit midpoint), a documented
second-order implicit Runge--Kutta method.
"""

from math import isfinite

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import (
    IntegrationConvergenceError,
    InvalidConfigurationError,
    NonPhysicalStateError,
)
from wildfire_model.models.parameters import ADRParameters
from wildfire_model.models.state import validate_temperature_biomass


def _active_midpoint_temperature(
    temperature: NDArray[np.float64],
    biomass: NDArray[np.float64],
    time_step: float,
    parameters: ADRParameters,
    *,
    maximum_iterations: int,
    relative_tolerance: float,
) -> tuple[NDArray[np.float64], NDArray[np.bool_]]:
    """Solve the active implicit midpoint stage by vectorized Newton iteration."""

    stage = temperature.copy()
    converged = np.zeros(temperature.shape, dtype=np.bool_)
    half_step = 0.5 * time_step
    for _ in range(maximum_iterations):
        rate = parameters.pre_exponential_factor * np.exp(
            -parameters.activation_temperature / stage
        )
        denominator = 1.0 + half_step * rate
        residual = (
            stage
            - temperature
            + half_step * parameters.cooling_rate * (stage - parameters.ambient_temperature)
            - half_step * parameters.heat_release_temperature * biomass * rate / denominator
        )
        rate_derivative = rate * parameters.activation_temperature / stage**2
        derivative = (
            1.0
            + half_step * parameters.cooling_rate
            - half_step
            * parameters.heat_release_temperature
            * biomass
            * rate_derivative
            / denominator**2
        )
        if np.any(~np.isfinite(derivative)) or np.any(np.abs(derivative) < 1.0e-14):
            raise IntegrationConvergenceError(
                "Implicit midpoint Newton derivative became singular or non-finite."
            )
        update = residual / derivative
        next_stage = stage - update
        if np.any(~np.isfinite(next_stage)) or np.any(next_stage <= 0.0):
            raise IntegrationConvergenceError(
                "Implicit midpoint Newton iteration proposed a nonphysical temperature."
            )
        converged_now = np.abs(update) <= relative_tolerance * (1.0 + np.abs(next_stage))
        converged |= converged_now
        stage = next_stage
        if np.all(converged):
            break
    return stage, converged


def implicit_midpoint_step(
    temperature: ArrayLike,
    biomass: ArrayLike,
    time_step: float,
    parameters: ADRParameters,
    *,
    maximum_iterations: int = 30,
    relative_tolerance: float = 1.0e-12,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Advance the local reaction system by one implicit midpoint step.

    The discontinuous ignition switch is evaluated at the Runge--Kutta
    midpoint. If the active-stage solution falls below ``T_pc``, the
    exact inactive midpoint stage (ambient cooling only) is used.
    """

    temperature_array, biomass_array = validate_temperature_biomass(
        temperature, biomass, context="implicit midpoint input"
    )
    if not isfinite(time_step) or time_step <= 0.0:
        raise InvalidConfigurationError("Reaction time_step must be finite and > 0 s.")
    if maximum_iterations < 1:
        raise InvalidConfigurationError("maximum_iterations must be at least one.")
    if not isfinite(relative_tolerance) or relative_tolerance <= 0.0:
        raise InvalidConfigurationError("relative_tolerance must be finite and > 0.")

    half_step = 0.5 * time_step
    inactive_stage = (
        temperature_array + half_step * parameters.cooling_rate * parameters.ambient_temperature
    ) / (1.0 + half_step * parameters.cooling_rate)
    active_candidates = temperature_array >= parameters.pyrolysis_temperature
    midpoint_temperature = inactive_stage.copy()
    midpoint_rate = np.zeros_like(temperature_array)

    if np.any(active_candidates):
        active_stage, converged = _active_midpoint_temperature(
            temperature_array[active_candidates],
            biomass_array[active_candidates],
            time_step,
            parameters,
            maximum_iterations=maximum_iterations,
            relative_tolerance=relative_tolerance,
        )
        if not np.all(converged):
            count = int(np.count_nonzero(~converged))
            raise IntegrationConvergenceError(
                f"Implicit midpoint Newton iteration failed in {count} state(s) after "
                f"{maximum_iterations} iterations."
            )
        stage_is_active = active_stage >= parameters.pyrolysis_temperature
        inactive_is_consistent = (
            inactive_stage[active_candidates] < parameters.pyrolysis_temperature
        )
        if np.any(~stage_is_active & ~inactive_is_consistent):
            raise IntegrationConvergenceError(
                "Implicit midpoint found neither an active nor an inactive "
                "switch-consistent stage; reduce the time step."
            )
        active_indices = np.flatnonzero(active_candidates)
        selected_indices = active_indices[stage_is_active]
        midpoint_temperature.flat[selected_indices] = active_stage[stage_is_active]
        midpoint_rate.flat[selected_indices] = parameters.pre_exponential_factor * np.exp(
            -parameters.activation_temperature / midpoint_temperature.flat[selected_indices]
        )

    positivity_factor = half_step * midpoint_rate
    if np.any(positivity_factor >= 1.0):
        maximum = float(np.max(positivity_factor))
        raise NonPhysicalStateError(
            "Implicit midpoint would produce negative biomass because "
            f"0.5*dt*Psi(T_mid)={maximum:.16g} is not below one; reduce the time step."
        )
    next_biomass = biomass_array * (1.0 - positivity_factor) / (1.0 + positivity_factor)
    next_temperature = 2.0 * midpoint_temperature - temperature_array
    return validate_temperature_biomass(
        next_temperature, next_biomass, context="implicit midpoint output"
    )


def implicit_midpoint_step_constant_reaction(
    temperature: ArrayLike,
    biomass: ArrayLike,
    time_step: float,
    parameters: ADRParameters,
    reaction_rate: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Advance Navas-Montilla Eqs. (41)-(42) by implicit midpoint."""

    temperature_array, biomass_array = validate_temperature_biomass(
        temperature, biomass, context="constant-reaction midpoint input"
    )
    if not isfinite(time_step) or time_step <= 0.0:
        raise InvalidConfigurationError("Reaction time_step must be finite and > 0 s.")
    if not isfinite(reaction_rate) or reaction_rate < 0.0:
        raise InvalidConfigurationError("reaction_rate must be finite and non-negative.")
    half_step = 0.5 * time_step
    active = temperature_array >= parameters.pyrolysis_temperature
    midpoint_rate = np.where(active, reaction_rate, 0.0)
    denominator = 1.0 + half_step * midpoint_rate
    midpoint_biomass = biomass_array / denominator
    midpoint_temperature = (
        temperature_array
        + half_step
        * (
            parameters.cooling_rate * parameters.ambient_temperature
            + parameters.heat_release_temperature * midpoint_rate * midpoint_biomass
        )
    ) / (1.0 + half_step * parameters.cooling_rate)
    active = active & (midpoint_temperature >= parameters.pyrolysis_temperature)
    midpoint_rate = np.where(active, reaction_rate, 0.0)
    denominator = 1.0 + half_step * midpoint_rate
    midpoint_biomass = biomass_array / denominator
    midpoint_temperature = (
        temperature_array
        + half_step
        * (
            parameters.cooling_rate * parameters.ambient_temperature
            + parameters.heat_release_temperature * midpoint_rate * midpoint_biomass
        )
    ) / (1.0 + half_step * parameters.cooling_rate)
    positivity_factor = half_step * midpoint_rate
    if np.any(positivity_factor >= 1.0):
        raise NonPhysicalStateError("Implicit midpoint constant reaction requires 0.5*dt*A < 1.")
    next_temperature = 2.0 * midpoint_temperature - temperature_array
    next_biomass = biomass_array * (1.0 - positivity_factor) / denominator
    return validate_temperature_biomass(
        next_temperature,
        next_biomass,
        context="constant-reaction midpoint output",
    )


def integrate_reaction(
    temperature: ArrayLike,
    biomass: ArrayLike,
    final_time: float,
    time_step: float,
    parameters: ADRParameters,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Integrate to ``final_time`` with deterministic final-step shortening."""

    if not isfinite(final_time) or final_time < 0.0:
        raise InvalidConfigurationError("final_time must be finite and >= 0 s.")
    if not isfinite(time_step) or time_step <= 0.0:
        raise InvalidConfigurationError("time_step must be finite and > 0 s.")
    current_temperature, current_biomass = validate_temperature_biomass(
        temperature, biomass, context="reaction integration input"
    )
    current_temperature = current_temperature.copy()
    current_biomass = current_biomass.copy()
    elapsed = 0.0
    while elapsed < final_time:
        step = min(time_step, final_time - elapsed)
        current_temperature, current_biomass = implicit_midpoint_step(
            current_temperature, current_biomass, step, parameters
        )
        elapsed += step
    return current_temperature, current_biomass
