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

The code uses finite volumes: it tracks average values in each cell and the
flux of heat through cell faces. WENO reconstructs values for advection,
central differences handle diffusion, and a CFL limit chooses a stable time
step. Every 2D field has shape `(ny, nx)`, with y on axis 0 and x on axis 1.
Inputs use SI units. [The technical model note](docs/technical_model.md)
records the equations and numerical choices, and
[open modeling questions](docs/open_modeling_questions.md) records source
ambiguities.

## Code layout

| Path | Purpose |
|---|---|
| `src/wildfire_model/models/` | Combustion, moisture, and material properties |
| `src/wildfire_model/numerics/` | Fluxes, time integration, boundaries, and solvers |
| `src/wildfire_model/mesh/` | Uniform 1D and 2D meshes |
| `src/wildfire_model/inputs/` | Configurations, fields, and initial conditions |
| `src/wildfire_model/helpers/` | Front and burned-area measurements |
| `src/wildfire_model/main.py` | Public simulation entry points |
| `examples/wildfire_model.ipynb` | Guided 2D simulation with inline results |

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

For the 1D Reisch case, use:

```bash
python scripts/run_adr_1d.py
```

This reads `src/wildfire_model/inputs/reisch_stage1.toml` and prints its
final values. Pass `--config path/to/case.toml` to run another 1D TOML
configuration. The file sets physical parameters, the mesh, numerical
settings, final time, and output interval.
