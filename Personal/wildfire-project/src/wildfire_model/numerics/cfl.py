"""Courant-Friedrichs-Lewy time-step control."""

from math import isfinite

import numpy as np
from numpy.typing import ArrayLike

from wildfire_model.exceptions import InvalidConfigurationError
from wildfire_model.mesh.mesh_1d import UniformMesh1D
from wildfire_model.mesh.mesh_2d import UniformMesh2D


def stable_transport_timestep(
    mesh: UniformMesh1D,
    velocity: float,
    diffusivity: float,
    cfl: float,
) -> float:
    r"""Return the combined explicit advection-diffusion stability bound.

    ``dt = CFL / (|v|/dx + 2 D/dx^2)``.

    The factor two is the one-dimensional forward-Euler central-
    diffusion bound. Reisch et al. state dynamic CFL control but do not
    publish a formula; this conservative additive spectral-rate bound is
    the documented Stage-1 decision. If both rates vanish, the returned
    value is infinity because transport imposes no restriction.
    """

    if not isfinite(velocity):
        raise InvalidConfigurationError("velocity must be finite.")
    if not isfinite(diffusivity) or diffusivity < 0.0:
        raise InvalidConfigurationError("diffusivity must be finite and >= 0.")
    if not isfinite(cfl) or not 0.0 < cfl <= 1.0:
        raise InvalidConfigurationError("cfl must lie in (0, 1].")
    rate = abs(velocity) / mesh.spacing + 2.0 * diffusivity / mesh.spacing**2
    return float("inf") if rate == 0.0 else cfl / rate


def stable_transport_timestep_2d(
    mesh: UniformMesh2D,
    velocity_x: ArrayLike,
    velocity_y: ArrayLike,
    diffusivity: ArrayLike,
    cfl: float,
) -> float:
    r"""Return the explicit 2D advection-diffusion stability bound.

    The project uses the conservative additive rate

    ``max(|u|/dx + |v|/dy + 2 D/dx² + 2 D/dy²)``

    because the referenced papers state dynamic CFL control without
    publishing a combined two-dimensional formula.
    """

    if not isfinite(cfl) or not 0.0 < cfl <= 1.0:
        raise InvalidConfigurationError("cfl must lie in (0, 1].")

    def checked_field(values: ArrayLike, name: str) -> np.ndarray:
        array = np.asarray(values, dtype=np.float64)
        if array.ndim == 0:
            array = np.full(mesh.shape, float(array), dtype=np.float64)
        if array.shape != mesh.shape or not np.all(np.isfinite(array)):
            raise InvalidConfigurationError(
                f"{name} must be finite and scalar or have shape {mesh.shape}."
            )
        return array

    x_velocity = checked_field(velocity_x, "x velocity")
    y_velocity = checked_field(velocity_y, "y velocity")
    coefficient = checked_field(diffusivity, "diffusivity")
    if np.any(coefficient < 0.0):
        raise InvalidConfigurationError("diffusivity must be non-negative.")
    local_rate = (
        np.abs(x_velocity) / mesh.dx
        + np.abs(y_velocity) / mesh.dy
        + 2.0 * coefficient / mesh.dx**2
        + 2.0 * coefficient / mesh.dy**2
    )
    maximum_rate = float(np.max(local_rate))
    return float("inf") if maximum_rate == 0.0 else cfl / maximum_rate
