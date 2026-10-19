"""Explicit time integration for the transport operator."""

from collections.abc import Callable
from math import isfinite

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError

ExplicitOperator = Callable[[NDArray[np.float64]], NDArray[np.float64]]
TimeDependentOperator = Callable[[float, NDArray[np.float64]], NDArray[np.float64]]


def ssprk3_step(
    values: ArrayLike,
    time_step: float,
    operator: ExplicitOperator,
) -> NDArray[np.float64]:
    """Advance one Shu-Osher SSP-RK3 step.

    This is the method listed by Reisch et al. Sec. 2.2 and written
    explicitly by Navas-Montilla et al. (2025), Sec. 3.2, Eq. (78).
    """

    initial = np.asarray(values, dtype=np.float64)
    if initial.ndim < 1 or initial.size == 0 or not np.all(np.isfinite(initial)):
        raise InvalidConfigurationError("SSP-RK3 input must be a finite non-empty array.")
    if not isfinite(time_step) or time_step <= 0.0:
        raise InvalidConfigurationError("SSP-RK3 time_step must be finite and > 0 s.")
    stage_one = initial + time_step * operator(initial)
    stage_two = 0.75 * initial + 0.25 * (stage_one + time_step * operator(stage_one))
    result = (1.0 / 3.0) * initial + (2.0 / 3.0) * (stage_two + time_step * operator(stage_two))
    if not np.all(np.isfinite(result)):
        raise InvalidConfigurationError("SSP-RK3 produced NaN or infinity.")
    return result


def ssprk3_step_time_dependent(
    values: ArrayLike,
    time: float,
    time_step: float,
    operator: TimeDependentOperator,
) -> NDArray[np.float64]:
    """Advance non-autonomous transport with third-order SSP-RK3."""

    initial = np.asarray(values, dtype=np.float64)
    if initial.ndim < 1 or initial.size == 0 or not np.all(np.isfinite(initial)):
        raise InvalidConfigurationError("SSP-RK3 input must be a finite non-empty array.")
    if not isfinite(time):
        raise InvalidConfigurationError("SSP-RK3 time must be finite.")
    if not isfinite(time_step) or time_step <= 0.0:
        raise InvalidConfigurationError("SSP-RK3 time_step must be finite and > 0 s.")
    stage_one = initial + time_step * operator(time, initial)
    stage_two = 0.75 * initial + 0.25 * (
        stage_one + time_step * operator(time + time_step, stage_one)
    )
    result = (1.0 / 3.0) * initial + (2.0 / 3.0) * (
        stage_two + time_step * operator(time + 0.5 * time_step, stage_two)
    )
    if not np.all(np.isfinite(result)):
        raise InvalidConfigurationError("SSP-RK3 produced NaN or infinity.")
    return result
