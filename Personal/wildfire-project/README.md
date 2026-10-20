# Wildfire model

This project simulates how heat and fuel change as a fire spreads. It has a
one-dimensional model based on Reisch et al. and a two-dimensional model
based on Navas-Montilla et al. Both are deterministic: the same inputs give
the same result.

## What the model computes

Each mesh cell stores temperature `T` in kelvin and the fraction of dry fuel
remaining `Y`, from 1 (unburned) to 0 (consumed). Heat moves to nearby cells
by diffusion and is carried by wind or slope-driven flow. When a cell reaches
its ignition temperature, fuel burns: `Y` falls and combustion adds heat.
Heat also escapes to the surroundings. This feedback can sustain a moving
fire front.

The 1D model uses fixed material properties. The complete 2D model can vary
fuel, moisture, wind, and terrain over space. It also changes the effective
heat capacity as water evaporates and adds a local radiative contribution to
heat transport. The supported moisture choices are `NM` (dry), `S2M`, and
`S3M`. `C2M` is held back because required values are missing from the
source paper.

Terrain is represented in 2.5D: the fire evolves over horizontal cells and
the result is drawn on a curved ground surface. Slope changes the effective
transport velocity. The [source model](https://arxiv.org/html/2412.04517v2)
averages over vegetation height, so its equations do not resolve vertical
flames or smoke.

The code uses finite volumes: it tracks average values in each cell and the
flux of heat through cell faces. WENO reconstructs values for advection,
central differences handle diffusion, and a CFL limit chooses a stable time
step. Every 2D field has shape `(ny, nx)`, with y on axis 0 and x on axis 1.
Inputs use SI units. The papers listed below give the equations and
numerical methods used by the project.

## Code layout

| Path | Purpose |
|---|---|
| `src/wildfire_model/models/` | Combustion, moisture, and material properties |
| `src/wildfire_model/numerics/` | Fluxes, time integration, boundaries, and solvers |
| `src/wildfire_model/mesh/` | Uniform 1D and 2D meshes |
| `src/wildfire_model/inputs/` | Configurations, fields, and initial conditions |
| `src/wildfire_model/helpers/` | Diagnostics and ParaView export |
| `src/wildfire_model/main.py` | Public simulation entry points |
| `examples/wildfire_model.ipynb` | Guided 2D simulation with inline results |
| `downloads/` | Generated XDMF/HDF5 results for ParaView |

The main outputs are temperature and remaining-fuel arrays at requested
times. A 2D result has shape `(time, ny, nx)` for each quantity.

## Run it

From a fresh clone, use Python 3.12 or newer:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[notebook]'
jupyter lab examples/wildfire_model.ipynb
```

Run the notebook cells from top to bottom. It sets up a small 2D fire,
advances the model, and shows temperature, remaining fuel, and a few
measurements directly in the notebook. It does not write figures to disk.
It also exports its terrain time series to `downloads/notebook/`.

To generate the terrain animation data directly:

```bash
python scripts/run_terrain_fire.py
```

This saves 121 snapshots over 180 seconds to
`downloads/terrain_fire/series.xdmf` and `series.h5`. Open the `.xdmf` file
in ParaView, click **Apply**, color by `temperature_K` or `remaining_fuel`,
and press **Play**. The terrain already has its actual elevation. Use
**File → Save Animation** to record the sequence; keep the `.h5` file next
to the `.xdmf` file. These are [ParaView's animation controls](https://docs.paraview.org/en/latest/UsersGuide/animation.html).

For the 1D Reisch case, use:

```bash
python scripts/run_adr_1d.py
```

This reads `src/wildfire_model/inputs/reisch_stage1.toml` and prints its
final values. Pass `--config path/to/case.toml` to run another 1D TOML
configuration. The file sets physical parameters, the mesh, numerical
settings, final time, and output interval.
The 1D command also exports `downloads/reisch_1d/series.xdmf` and `series.h5`.
Public simulation entry points export automatically. For a solver used
directly, `write_paraview_series` exports its result with a chosen terrain
and a directory under `downloads/`. Pass `output_directory=None` to a public
run function for an in-memory run.

## References

- Reisch, Navas-Montilla, and Özgen-Xian (2024).
  [Analytical and numerical insights into wildfire dynamics: Exploring the
  advection–diffusion–reaction model](https://doi.org/10.1016/j.camwa.2024.01.024).
  The dry ADR equations and 1D reference cases.
- Navas-Montilla, Reisch, Diaz, and Özgen-Xian (2025).
  [Modeling wildfire dynamics through a physics-based approach incorporating
  fuel moisture and landscape heterogeneity](https://arxiv.org/abs/2412.04517v2).
  The mixture, moisture, radiation, wind, and slope terms.
- Asensio, Cascón, Laiz, and Prieto-Herráez (2023).
  [Validating the effect of fuel moisture content by a multivalued operator
  in a simplified physical fire spread model](https://doi.org/10.1016/j.envsoft.2023.105710).
  Background on moisture and spread trends.

Output follows the [XDMF format specification](https://www.xdmf.org/index.php/XDMF_Model_and_Format.html)
using the [h5py dataset API](https://docs.h5py.org/en/stable/high/dataset.html).
