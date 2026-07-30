"""Run the Reisch wave spatial-resolution and propagation-speed study."""

from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from wildfire_model.helpers.postprocessing import front_position, propagation_speed
from wildfire_model.inputs.reference_cases import reisch_wave_configuration
from wildfire_model.main import run_simulation


def main() -> None:
    output_directory = Path("artifacts/convergence_study")
    output_directory.mkdir(parents=True, exist_ok=True)
    cell_counts = (400, 800, 1600)
    speeds = []
    for n_cells in cell_counts:
        configuration = replace(
            reisch_wave_configuration(n_cells),
            output_times=np.array([0.0, 50.0, 100.0, 150.0, 200.0]),
        )
        simulation = run_simulation(configuration)
        result = simulation.result
        positions = np.array(
            [
                front_position(
                    simulation.mesh.cell_centers,
                    temperature,
                    threshold=400.0,
                    direction="right",
                )
                for temperature in result.temperature
            ]
        )
        speeds.append(propagation_speed(result.times[1:], positions[1:]))

    spacings = 500.0 / np.asarray(cell_counts)
    with (output_directory / "wave_speed_convergence.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(["cells", "dx_m", "right_front_speed_m_per_s"])
        writer.writerows(zip(cell_counts, spacings, speeds, strict=True))

    figure, axis = plt.subplots(figsize=(7, 4.5))
    axis.plot(spacings, speeds, "o-")
    axis.invert_xaxis()
    axis.set_xlabel("Cell width (m)")
    axis.set_ylabel("Right-front speed (m/s)")
    axis.grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(output_directory / "wave_speed_convergence.png", dpi=160)
    plt.close(figure)

    for n_cells, speed in zip(cell_counts, speeds, strict=True):
        print(f"cells={n_cells} right_front_speed_m_per_s={speed:.12g}")


if __name__ == "__main__":
    main()
