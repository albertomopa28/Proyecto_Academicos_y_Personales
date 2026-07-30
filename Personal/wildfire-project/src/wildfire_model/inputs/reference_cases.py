"""Published Reisch verification and wave cases in strict SI units."""

import numpy as np
from numpy.typing import NDArray
from scipy.special import erf

from wildfire_model.inputs.config import MeshConfiguration, SimulationConfiguration
from wildfire_model.inputs.parameters import (
    reisch_case_2_parameters,
    reisch_reaction_parameters,
)
from wildfire_model.mesh.mesh_1d import UniformMesh1D
from wildfire_model.models.state import ADRState


def _gaussian_cell_averages(
    mesh: UniformMesh1D,
    *,
    centre: float,
    coefficient: float,
    amplitude: float,
    baseline: float,
) -> NDArray[np.float64]:
    faces = mesh.faces
    scale = np.sqrt(coefficient)
    integral = (
        0.5
        * np.sqrt(np.pi / coefficient)
        * (erf(scale * (faces[1:] - centre)) - erf(scale * (faces[:-1] - centre)))
    )
    return np.asarray(baseline + amplitude * integral / mesh.spacing, dtype=np.float64)


def reisch_case_2_exact_cell_averages(
    mesh: UniformMesh1D,
    time: float,
) -> NDArray[np.float64]:
    """Cell averages of Reisch Appendix A.2, Eqs. (44)--(45)."""

    parameters = reisch_case_2_parameters()
    denominator = 1.0 + 0.004 * parameters.thermal_diffusivity * time
    return _gaussian_cell_averages(
        mesh,
        centre=250.0 + parameters.velocity * time,
        coefficient=0.001 / denominator,
        amplitude=100.0 * np.exp(-parameters.cooling_rate * time) / np.sqrt(denominator),
        baseline=parameters.ambient_temperature,
    )


def reisch_case_2_initial_state(mesh: UniformMesh1D) -> ADRState:
    """Initial cell averages for Reisch Appendix A.2."""

    return ADRState.from_arrays(
        reisch_case_2_exact_cell_averages(mesh, 0.0),
        np.ones(mesh.n_cells, dtype=np.float64),
    )


def reisch_wave_initial_state(mesh: UniformMesh1D) -> ADRState:
    """Reisch Sec. 3.3 Eq. (36), represented as exact cell averages."""

    faces = mesh.faces
    overlap = np.maximum(
        0.0,
        np.minimum(faces[1:], 275.0) - np.maximum(faces[:-1], 225.0),
    )
    hot_fraction = overlap / mesh.spacing
    temperature = 300.0 + 170.0 * hot_fraction
    return ADRState.from_arrays(
        temperature,
        np.ones(mesh.n_cells, dtype=np.float64),
    )


def reisch_wave_configuration(n_cells: int = 800) -> SimulationConfiguration:
    """Return the Reisch Eq. (36) wave configuration used for verification."""

    return SimulationConfiguration(
        parameters=reisch_reaction_parameters(),
        mesh=MeshConfiguration(0.0, 500.0, n_cells),
        cfl=0.1,
        weno_order=7,
        boundary="transmissive",
        output_times=np.linspace(0.0, 200.0, 9),
        maximum_time_step=0.05,
        initial_condition="reisch_wave",
    )
