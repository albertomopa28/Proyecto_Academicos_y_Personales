r"""Finite-volume WENO-JS reconstruction of orders 1, 3, 5, and 7.

Reisch et al. (2024), Sec. 2.2, specify finite volumes with WENO up to
seventh order and cite Jiang & Shu (1996). The implementation below
constructs the Jiang-Shu candidate polynomials directly from cell
averages. Smoothness indicators integrate squared derivatives over the
upwind cell. This avoids separate, error-prone hard-coded formulae for
each supported order.
"""

from dataclasses import dataclass
from functools import lru_cache
from math import isfinite

import numpy as np
from numpy.typing import ArrayLike, NDArray

from wildfire_model.exceptions import InvalidConfigurationError
from wildfire_model.numerics.boundary_conditions import (
    BoundaryCondition,
    extend_with_ghost_cells,
)

SUPPORTED_WENO_ORDERS = (1, 3, 5, 7)


@dataclass(frozen=True, slots=True)
class _WENOData:
    radius: int
    offsets: NDArray[np.int64]
    candidate_indices: tuple[NDArray[np.int64], ...]
    candidate_coefficients: tuple[NDArray[np.float64], ...]
    smoothness_matrices: tuple[NDArray[np.float64], ...]
    linear_weights: NDArray[np.float64]


def _average_matrix(offsets: NDArray[np.int64]) -> NDArray[np.float64]:
    degree = offsets.size - 1
    matrix = np.empty((offsets.size, degree + 1), dtype=np.float64)
    for row, offset in enumerate(offsets):
        left = float(offset) - 0.5
        right = float(offset) + 0.5
        for power in range(degree + 1):
            matrix[row, power] = (right ** (power + 1) - left ** (power + 1)) / (power + 1)
    return matrix


def _evaluation_coefficients(offsets: NDArray[np.int64]) -> NDArray[np.float64]:
    matrix = _average_matrix(offsets)
    powers = np.power(0.5, np.arange(offsets.size, dtype=np.float64))
    return np.linalg.solve(matrix.T, powers)


def _smoothness_matrix(offsets: NDArray[np.int64]) -> NDArray[np.float64]:
    radius = offsets.size
    if radius == 1:
        return np.zeros((1, 1), dtype=np.float64)
    polynomial_matrix = np.linalg.inv(_average_matrix(offsets))
    derivative_energy = np.zeros((radius, radius), dtype=np.float64)
    for derivative_order in range(1, radius):
        for first_power in range(derivative_order, radius):
            first_factor = float(
                np.prod(np.arange(first_power - derivative_order + 1, first_power + 1))
            )
            for second_power in range(derivative_order, radius):
                second_factor = float(
                    np.prod(
                        np.arange(
                            second_power - derivative_order + 1,
                            second_power + 1,
                        )
                    )
                )
                exponent = first_power + second_power - 2 * derivative_order
                integral = (0.5 ** (exponent + 1) - (-0.5) ** (exponent + 1)) / (exponent + 1)
                derivative_energy[first_power, second_power] += (
                    first_factor * second_factor * integral
                )
    return polynomial_matrix.T @ derivative_energy @ polynomial_matrix


@lru_cache(maxsize=4)
def _weno_data(order: int) -> _WENOData:
    if order not in SUPPORTED_WENO_ORDERS:
        raise InvalidConfigurationError(
            f"WENO order must be one of {SUPPORTED_WENO_ORDERS}; received {order}."
        )
    radius = (order + 1) // 2
    offsets = np.arange(-radius + 1, radius, dtype=np.int64)
    candidate_indices: list[NDArray[np.int64]] = []
    candidate_coefficients: list[NDArray[np.float64]] = []
    smoothness_matrices: list[NDArray[np.float64]] = []
    embedded_coefficients: list[NDArray[np.float64]] = []
    for candidate in range(radius):
        indices = np.arange(candidate, candidate + radius, dtype=np.int64)
        local_offsets = offsets[indices]
        coefficients = _evaluation_coefficients(local_offsets)
        embedded = np.zeros(order, dtype=np.float64)
        embedded[indices] = coefficients
        candidate_indices.append(indices)
        candidate_coefficients.append(coefficients)
        smoothness_matrices.append(_smoothness_matrix(local_offsets))
        embedded_coefficients.append(embedded)
    full_coefficients = _evaluation_coefficients(offsets)
    candidate_matrix = np.stack(embedded_coefficients, axis=1)
    linear_weights, residuals, rank, _ = np.linalg.lstsq(
        candidate_matrix, full_coefficients, rcond=None
    )
    if (
        rank != radius
        or (residuals.size and float(residuals[0]) > 1.0e-24)
        or np.any(linear_weights <= 0.0)
    ):
        raise InvalidConfigurationError(
            f"Failed to construct positive optimal WENO-{order} weights."
        )
    return _WENOData(
        radius=radius,
        offsets=offsets,
        candidate_indices=tuple(candidate_indices),
        candidate_coefficients=tuple(candidate_coefficients),
        smoothness_matrices=tuple(smoothness_matrices),
        linear_weights=linear_weights,
    )


def _reconstruct_stencil(
    stencil: NDArray[np.float64],
    data: _WENOData,
    epsilon: float,
) -> float:
    if data.radius == 1:
        return float(stencil[0])
    candidate_values = np.empty(data.radius, dtype=np.float64)
    smoothness = np.empty(data.radius, dtype=np.float64)
    for index in range(data.radius):
        local_values = stencil[data.candidate_indices[index]]
        candidate_values[index] = data.candidate_coefficients[index] @ local_values
        smoothness[index] = local_values @ data.smoothness_matrices[index] @ local_values
    nonlinear_weights = data.linear_weights / (epsilon + smoothness) ** 2
    nonlinear_weights /= np.sum(nonlinear_weights)
    return float(nonlinear_weights @ candidate_values)


def reconstruct_faces(
    values: ArrayLike,
    order: int,
    boundary: BoundaryCondition,
    *,
    epsilon: float = 1.0e-6,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Reconstruct left and right states at all ``N+1`` cell faces."""

    if not isfinite(epsilon) or epsilon <= 0.0:
        raise InvalidConfigurationError("WENO epsilon must be finite and > 0.")
    array = np.asarray(values, dtype=np.float64)
    if array.ndim != 1 or array.size == 0 or not np.all(np.isfinite(array)):
        raise InvalidConfigurationError("WENO input must be a finite non-empty 1D array.")
    data = _weno_data(order)
    extended = extend_with_ghost_cells(array, data.radius, boundary)
    left_states = np.empty(array.size + 1, dtype=np.float64)
    right_states = np.empty(array.size + 1, dtype=np.float64)
    for face in range(array.size + 1):
        left_anchor = face - 1 + data.radius
        left_stencil = extended[left_anchor + data.offsets]
        left_states[face] = _reconstruct_stencil(left_stencil, data, epsilon)

        right_anchor = face + data.radius
        right_stencil = extended[right_anchor - data.offsets]
        right_states[face] = _reconstruct_stencil(right_stencil, data, epsilon)
    return left_states, right_states


def upwind_flux(
    values: ArrayLike,
    velocity: float,
    order: int,
    boundary: BoundaryCondition,
) -> NDArray[np.float64]:
    """Return ``v*T`` at faces using the state upwind of each face."""

    if not isfinite(velocity):
        raise InvalidConfigurationError("Advection velocity must be finite.")
    array = np.asarray(values, dtype=np.float64)
    if velocity == 0.0:
        return np.zeros(array.size + 1, dtype=np.float64)
    left_states, right_states = reconstruct_faces(array, order, boundary)
    return velocity * (left_states if velocity > 0.0 else right_states)
