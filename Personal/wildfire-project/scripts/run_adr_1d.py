"""Run a deterministic one-dimensional ADR simulation from TOML input."""

from __future__ import annotations

import argparse
from pathlib import Path

from wildfire_model.main import run_from_config


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("src/wildfire_model/inputs/reisch_stage1.toml"),
        help="TOML simulation configuration.",
    )
    arguments = parser.parse_args()
    simulation = run_from_config(arguments.config)
    result = simulation.result
    print(f"cells={simulation.mesh.n_cells}")
    print(f"final_time_s={result.times[-1]:.12g}")
    print(f"maximum_temperature_K={result.temperature[-1].max():.12g}")
    print(f"minimum_biomass={result.biomass[-1].min():.12g}")


if __name__ == "__main__":
    main()
