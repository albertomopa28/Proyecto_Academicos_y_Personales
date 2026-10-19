"""Finite-volume advection and second-order diffusion operators."""

from math import isfinite

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError
from wildfire_model.mesh.mesh_1d import UniformMesh1D
from wildfire_model.mesh.mesh_2d import UniformMesh2D
from wildfire_model.numerics.boundary_conditions import (
    BoundaryCondition,
    BoundaryCondition2D,
    extend_with_ghost_cells,
    extend_with_ghost_cells_2d,
    reconstruction_boundary,
)
from wildfire_model.numerics.weno import (
    reconstruct_faces_along_axis,
    upwind_flux,
)


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


def _validate_2d_field(
    values: ArrayLike,
    mesh: UniformMesh2D,
    name: str,
    *,
    nonnegative: bool = False,
) -> NDArray[np.float64]:
    array = np.asarray(values, dtype=np.float64)
    if array.ndim == 0:
        array = np.full(mesh.shape, float(array), dtype=np.float64)
    if array.shape != mesh.shape or not np.all(np.isfinite(array)):
        raise InvalidConfigurationError(
            f"{name} must be finite with shape {mesh.shape}; received {array.shape}."
        )
    if nonnegative and np.any(array < 0.0):
        raise InvalidConfigurationError(f"{name} must be non-negative.")
    return array


def _face_averages_2d(
    values: NDArray[np.float64],
    boundary: BoundaryCondition2D,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    extended = extend_with_ghost_cells_2d(values, 1, boundary)
    x_faces = 0.5 * (extended[1:-1, :-1] + extended[1:-1, 1:])
    y_faces = 0.5 * (extended[:-1, 1:-1] + extended[1:, 1:-1])
    return x_faces, y_faces


def advection_fluxes_2d(
    values: ArrayLike,
    mesh: UniformMesh2D,
    velocity_x: ArrayLike,
    velocity_y: ArrayLike,
    order: int,
    boundary: BoundaryCondition2D,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return conservative upwind ``v*T`` fluxes on all Cartesian faces."""

    array = _validate_2d_field(values, mesh, "Advection input")
    velocity_x_array = _validate_2d_field(velocity_x, mesh, "x velocity")
    velocity_y_array = _validate_2d_field(velocity_y, mesh, "y velocity")
    ghost_boundary = reconstruction_boundary(boundary)
    left, right = reconstruct_faces_along_axis(array, order, ghost_boundary, axis=1)
    lower, upper = reconstruct_faces_along_axis(array, order, ghost_boundary, axis=0)
    face_velocity_x, face_velocity_y = (
        _face_averages_2d(velocity_x_array, boundary)[0],
        _face_averages_2d(velocity_y_array, boundary)[1],
    )
    flux_x = face_velocity_x * np.where(face_velocity_x >= 0.0, left, right)
    flux_y = face_velocity_y * np.where(face_velocity_y >= 0.0, lower, upper)
    return flux_x, flux_y


def diffusive_fluxes_2d(
    values: ArrayLike,
    mesh: UniformMesh2D,
    diffusivity: ArrayLike,
    boundary: BoundaryCondition2D,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return ``D grad(T)`` using Navas-Montilla Eq. (76) face means."""

    array = _validate_2d_field(values, mesh, "Diffusion input")
    coefficient = _validate_2d_field(diffusivity, mesh, "Thermal diffusivity", nonnegative=True)
    extended = extend_with_ghost_cells_2d(array, 1, boundary)
    coefficient_x, coefficient_y = _face_averages_2d(coefficient, boundary)
    gradient_x = (extended[1:-1, 1:] - extended[1:-1, :-1]) / mesh.dx
    gradient_y = (extended[1:, 1:-1] - extended[:-1, 1:-1]) / mesh.dy
    flux_x = coefficient_x * gradient_x
    flux_y = coefficient_y * gradient_y
    if boundary != "periodic":
        flux_x[:, (0, -1)] = 0.0
        flux_y[(0, -1), :] = 0.0
    return flux_x, flux_y


def transport_rhs_2d(
    values: ArrayLike,
    mesh: UniformMesh2D,
    velocity_x: ArrayLike,
    velocity_y: ArrayLike,
    diffusivity: ArrayLike,
    order: int,
    boundary: BoundaryCondition2D,
) -> NDArray[np.float64]:
    """Return the conservative two-dimensional advection-diffusion operator."""

    advective_x, advective_y = advection_fluxes_2d(
        values, mesh, velocity_x, velocity_y, order, boundary
    )
    diffusive_x, diffusive_y = diffusive_fluxes_2d(values, mesh, diffusivity, boundary)
    total_x = advective_x - diffusive_x
    total_y = advective_y - diffusive_y
    if boundary == "zero_total_flux":
        total_x[:, (0, -1)] = 0.0
        total_y[(0, -1), :] = 0.0
    return -(
        (total_x[:, 1:] - total_x[:, :-1]) / mesh.dx + (total_y[1:, :] - total_y[:-1, :]) / mesh.dy
    )


def energy_transport_rhs_2d(
    temperature: ArrayLike,
    mesh: UniformMesh2D,
    velocity_x: ArrayLike,
    velocity_y: ArrayLike,
    thermal_conductivity: ArrayLike,
    volumetric_heat_capacity: ArrayLike,
    order: int,
    boundary: BoundaryCondition2D,
) -> NDArray[np.float64]:
    """Discretize the transport terms of Navas-Montilla Eq. (32).

    Conductive and radiative fluxes use Eq. (76), then their divergence is
    divided by cell-centered volumetric heat capacity. Advection uses the
    conservative flux form selected for the two-dimensional core.
    """

    capacity = _validate_2d_field(
        volumetric_heat_capacity,
        mesh,
        "Volumetric heat capacity",
    )
    if np.any(capacity <= 0.0):
        raise InvalidConfigurationError("Volumetric heat capacity must be strictly positive.")
    advective_x, advective_y = advection_fluxes_2d(
        temperature,
        mesh,
        velocity_x,
        velocity_y,
        order,
        boundary,
    )
    conductive_x, conductive_y = diffusive_fluxes_2d(
        temperature,
        mesh,
        thermal_conductivity,
        boundary,
    )
    if boundary == "zero_total_flux":
        capacity_x, capacity_y = _face_averages_2d(capacity, boundary)
        total_x = advective_x - conductive_x / capacity_x
        total_y = advective_y - conductive_y / capacity_y
        total_x[:, (0, -1)] = 0.0
        total_y[(0, -1), :] = 0.0
        return -(
            (total_x[:, 1:] - total_x[:, :-1]) / mesh.dx
            + (total_y[1:, :] - total_y[:-1, :]) / mesh.dy
        )
    advection = -(
        (advective_x[:, 1:] - advective_x[:, :-1]) / mesh.dx
        + (advective_y[1:, :] - advective_y[:-1, :]) / mesh.dy
    )
    conduction = (
        (conductive_x[:, 1:] - conductive_x[:, :-1]) / mesh.dx
        + (conductive_y[1:, :] - conductive_y[:-1, :]) / mesh.dy
    ) / capacity
    return advection + conduction
