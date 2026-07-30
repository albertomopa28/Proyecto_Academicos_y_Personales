"""Reproduce the Reisch wave and WENO verification figures."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from wildfire_model.helpers.postprocessing import (
    front_position,
    observed_orders,
    propagation_speed,
)
from wildfire_model.inputs.reference_cases import reisch_wave_configuration
from wildfire_model.main import run_simulation
from wildfire_model.numerics.weno import reconstruct_faces


def _smooth_periodic_cell_averages(n_cells: int) -> tuple[np.ndarray, np.ndarray]:
    faces = np.linspace(0.0, 2.0 * np.pi, n_cells + 1)
    spacing = faces[1] - faces[0]
    averages = (
        np.cos(faces[:-1])
        - np.cos(faces[1:])
        + 0.1 * (np.cos(2.0 * faces[:-1]) - np.cos(2.0 * faces[1:]))
    ) / spacing
    exact_faces = np.sin(faces) + 0.2 * np.sin(2.0 * faces)
    return averages, exact_faces


def _run_wave(output_directory: Path) -> float:
    simulation = run_simulation(reisch_wave_configuration())
    mesh = simulation.mesh
    result = simulation.result
    times = result.times
    positions = np.array(
        [
            front_position(
                mesh.cell_centers,
                temperature,
                threshold=400.0,
                direction="right",
            )
            for temperature in result.temperature
        ]
    )
    fit_mask = times >= 50.0
    speed = propagation_speed(times[fit_mask], positions[fit_mask])

    with (output_directory / "wave_front.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time_s", "right_front_m", "maximum_temperature_K", "minimum_biomass"])
        writer.writerows(
            zip(
                times,
                positions,
                np.max(result.temperature, axis=1),
                np.min(result.biomass, axis=1),
                strict=True,
            )
        )

    selected = [0, 2, 4, 6, 8]
    figure, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for index in selected:
        label = f"{times[index]:.0f} s"
        axes[0].plot(mesh.cell_centers, result.temperature[index], label=label)
        axes[1].plot(mesh.cell_centers, result.biomass[index], label=label)
    axes[0].axhline(400.0, color="black", linestyle=":", linewidth=1)
    axes[0].set_ylabel("Temperature (K)")
    axes[1].set_ylabel("Remaining biomass (-)")
    axes[1].set_xlabel("x (m)")
    axes[0].legend(ncol=3)
    axes[1].set_ylim(-0.02, 1.02)
    figure.suptitle("Reisch Eq. (36): WENO-7 Stage-1 evolution")
    figure.tight_layout()
    figure.savefig(output_directory / "temperature_biomass_evolution.png", dpi=160)
    plt.close(figure)

    line = np.polyfit(times[fit_mask], positions[fit_mask], 1)
    figure, axis = plt.subplots(figsize=(7, 4.5))
    axis.plot(times, positions, "o-", label="interpolated 400 K front")
    axis.plot(
        times[fit_mask],
        np.polyval(line, times[fit_mask]),
        "--",
        label=f"fit after 50 s: {speed:.5f} m/s",
    )
    axis.set_xlabel("Time (s)")
    axis.set_ylabel("Right-front position (m)")
    axis.legend()
    axis.grid(alpha=0.25)
    figure.tight_layout()
    figure.savefig(output_directory / "front_position_speed.png", dpi=160)
    plt.close(figure)
    return speed


def _run_weno_convergence(output_directory: Path) -> dict[int, float]:
    orders = (1, 3, 5, 7)
    cell_counts = np.array([40, 80, 160, 320])
    rows: list[tuple[int, int, float, float, float]] = []
    finest_orders: dict[int, float] = {}
    figure, axis = plt.subplots(figsize=(7, 5))
    for order in orders:
        errors = []
        for n_cells in cell_counts:
            averages, exact_faces = _smooth_periodic_cell_averages(int(n_cells))
            reconstructed, _ = reconstruct_faces(averages, order, "periodic")
            errors.append(float(np.mean(np.abs(reconstructed[:-1] - exact_faces[:-1]))))
        spacings = 2.0 * np.pi / cell_counts
        measured = observed_orders(spacings, errors)
        finest_orders[order] = float(measured[-1])
        row_orders = np.concatenate(([np.nan], measured))
        rows.extend(
            (
                order,
                int(n_cells),
                float(spacing),
                error,
                float(measured_order),
            )
            for n_cells, spacing, error, measured_order in zip(
                cell_counts,
                spacings,
                errors,
                row_orders,
                strict=True,
            )
        )
        axis.loglog(
            spacings,
            errors,
            "o-",
            label=f"WENO-{order}, finest p={measured[-1]:.3f}",
        )

    with (output_directory / "weno_convergence.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.writer(stream)
        writer.writerow(["weno_order", "cells", "dx", "l1_error", "observed_order"])
        writer.writerows(rows)

    axis.set_xlabel("Cell width")
    axis.set_ylabel("Interface reconstruction L1 error")
    axis.invert_xaxis()
    axis.grid(which="both", alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_directory / "weno_order_comparison.png", dpi=160)
    plt.close(figure)
    return finest_orders


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("artifacts/reisch_verification"),
        help="Directory for CSV diagnostics and PNG figures.",
    )
    arguments = parser.parse_args()
    arguments.output.mkdir(parents=True, exist_ok=True)

    speed = _run_wave(arguments.output)
    orders = _run_weno_convergence(arguments.output)
    summary = [
        f"right_front_speed_m_per_s={speed:.12g}",
        *(f"weno_{order}_finest_observed_order={value:.12g}" for order, value in orders.items()),
    ]
    (arguments.output / "summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))
    print(f"outputs={arguments.output.resolve()}")


if __name__ == "__main__":
    main()
