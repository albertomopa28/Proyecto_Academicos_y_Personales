"""Courant-Friedrichs-Lewy time-step control."""

from math import isfinite

from wildfire_model.exceptions import InvalidConfigurationError
from wildfire_model.mesh.mesh_1d import UniformMesh1D


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
