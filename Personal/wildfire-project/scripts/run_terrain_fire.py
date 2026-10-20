"""Run the terrain fire and write a ParaView animation time series."""

import argparse
from pathlib import Path

from wildfire_model.inputs.reference_cases_2d import terrain_fire_case
from wildfire_model.main import run_complete_case


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=80, help="Cells along x.")
    parser.add_argument("--ny", type=int, default=64, help="Cells along y.")
    parser.add_argument("--duration", type=float, default=180.0, help="Simulation time in s.")
    parser.add_argument("--frames", type=int, default=121, help="Number of saved time snapshots.")
    parser.add_argument("--output", type=Path, default=Path("downloads/terrain_fire"))
    arguments = parser.parse_args()
    case = terrain_fire_case(
        nx=arguments.nx, ny=arguments.ny, duration=arguments.duration, frames=arguments.frames
    )
    result = run_complete_case(case, output_directory=arguments.output)
    print(f"paraview_file={(arguments.output / 'series.xdmf').resolve()}")
    print(f"frames={result.times.size}; final_time_s={result.times[-1]:.12g}")
    print(f"maximum_temperature_K={result.temperature[-1].max():.12g}")
    print(f"minimum_biomass={result.biomass[-1].min():.12g}")


if __name__ == "__main__":
    main()
