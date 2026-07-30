"""Finite-volume advection and second-order diffusion operators."""

from math import isfinite

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError
from wildfire_model.mesh.mesh_1d import UniformMesh1D
from wildfire_model.numerics.boundary_conditions import (
    BoundaryCondition,
    extend_with_ghost_cells,
)
from wildfire_model.numerics.weno import upwind_flux


def advection_rhs(
    values: ArrayLike,
    mesh: UniformMesh1D,
    velocity: float,
    order: int,
    boundary: BoundaryCondition,
) -> NDArray[np.float64]:
    """Return the finite-volume approximation of ``-v dT/dx``."""

    flux = upwind_flux(values, velocity, order, boundary)
    return -(flux[1:] - flux[:-1]) / mesh.spacing


def diffusion_rhs(
    values: ArrayLike,
    mesh: UniformMesh1D,
    diffusivity: float,
    boundary: BoundaryCondition,
) -> NDArray[np.float64]:
    """Return the central approximation of ``D d²T/dx²``.

    This is the one-dimensional form of Reisch et al. Eq. (14), using
    their Sec. 2.2 second-order central diffusive flux and
    ``D = k/(rho_0 C)``.
    """

    if not isfinite(diffusivity) or diffusivity < 0.0:
        raise InvalidConfigurationError("Diffusivity must be finite and >= 0 m^2/s.")
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size != mesh.n_cells or not np.all(np.isfinite(array)):
        raise InvalidConfigurationError(
            "Diffusion input must be finite with one value per mesh cell."
        )
    if diffusivity == 0.0:
        return np.zeros_like(array)
    extended = extend_with_ghost_cells(array, 1, boundary)
    face_gradient = (extended[1:] - extended[:-1]) / mesh.spacing
    return diffusivity * (face_gradient[1:] - face_gradient[:-1]) / mesh.spacing


def transport_rhs(
    values: ArrayLike,
    mesh: UniformMesh1D,
    velocity: float,
    diffusivity: float,
    order: int,
    boundary: BoundaryCondition,
) -> NDArray[np.float64]:
    """Return the combined advection-diffusion semi-discrete operator."""

    return advection_rhs(values, mesh, velocity, order, boundary) + diffusion_rhs(
        values, mesh, diffusivity, boundary
    )
