# Wildfire ADR model

This repository implements a deterministic one-dimensional
advection-diffusion-reaction wildfire model based on Reisch,
Navas-Montilla, and Özgen-Xian (2024). The solver evolves temperature and
remaining biomass with finite volumes, WENO reconstruction, central
diffusion, SSP-RK3 transport integration, implicit midpoint reaction
integration, Strang splitting, and CFL-controlled time steps.

The current scope is the dry, deterministic, uniform Cartesian
one-dimensional model. Two-dimensional propagation, moisture, topography,
heterogeneous fuel, geospatial inputs, PINNs, inference, and uncertainty
quantification are not implemented.

All library quantities use SI units: m, s, kg, K, J, and W.

## Repository structure

```text
src/wildfire_model/
├── main.py
├── models/       # Physical parameters, state, combustion, and reaction equations
├── numerics/     # Finite volumes, WENO, integrators, CFL, boundaries, and solver
├── mesh/         # Uniform one-dimensional mesh geometry
├── inputs/       # Configuration and Reisch reference inputs
└── helpers/      # Front, speed, error, and convergence diagnostics
scripts/
├── run_adr_1d.py
├── run_reisch_verification.py
└── run_convergence_study.py
docs/
references/
```

## Setup

```bash
python -m pip install -e '.[dev]'
```

## Scripts

```bash
python scripts/run_adr_1d.py
python scripts/run_reisch_verification.py
python scripts/run_convergence_study.py
```

Generated figures and tabular diagnostics are written under `artifacts/`.
The physical formulation and numerical scheme are described in
`docs/model_specification.md` and `docs/numerical_method.md`.
