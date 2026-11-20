# 2D RANS solver

This folder is a small cell-centered finite-volume SIMPLE solver. Its first
run is a straight 20 mm by 4 mm rectangular channel, with no side outlet. It
uses first-order upwind convection, Rhie-Chow face fluxes, pressure correction,
and a two-equation k–omega SST closure. Run it from the workspace root:

```powershell
.\.venv\Scripts\python.exe rans_solver\solve.py
```

Edit `case.py` to change the mesh or inputs. Defaults use water properties at
20 C, 5% inlet turbulence intensity, and a 0.56 mm turbulence length scale.
These are starting assumptions because the project brief does not give water
temperature or inlet turbulence measurements. At 1 mm spacing the mesh is a
coarse screening mesh; use a wall-resolved mesh before treating SST wall
stresses as verified.

## What the code does

- `case.py` holds geometry, fluid, inlet, outlet, turbulence, and iteration
  inputs in SI units.
- `solve.py` builds internal and boundary faces, assembles upwind transport
  equations, uses Rhie-Chow corrected face fluxes in a SIMPLE pressure
  correction, and updates k, omega, and the SST eddy viscosity.
- Outputs are written to `rans_solver/results/` as a compressed NumPy file and
  a velocity/pressure plot. Console output includes residual and mass balance.

Run the fixed-geometry manifold reference with `--case manifold`:

```powershell
.\.venv\Scripts\python.exe rans_solver\solve.py --case manifold
```

That case uses the explicit wall faces in `mesh.py` and is the reference for
later Brinkman models. It saves separately under
`rans_solver/results/manifold_reference/`; the rectangle remains under
`rans_solver/results/rectangle_baseline/`.

## Brinkman and topology models

Use `--model A`, `--model B`, or `--model C` for density-defined Brinkman
manifolds. A uses a Brinkman momentum sink; B adds a linearly reconstructed
fluid/solid interface; C uses guarded quadratic reconstruction. For example:

```powershell
.\.venv\Scripts\python.exe rans_solver\solve.py --case manifold --model B
.\.venv\Scripts\python.exe -m rans_solver.compare_models
.\.venv\Scripts\python.exe -m rans_solver.optimize --steps 5 --model B
.\.venv\Scripts\python.exe -m rans_solver.plot_optimized
```

Model details, formulas, and known limits are in `BRINKMAN.md`. The comparison
script writes flow split and residuals to CSV. The optimizer uses coarse,
finite-difference sensitivities, fixed inlet/outlet and channel paths, a
connectivity check, warm-started solves, and a fixed fluid-volume constraint.
It uses an exploratory field-change tolerance; see `BRINKMAN.md` before
interpreting its output. No 3D verification is included.
`plot_optimized` reruns the saved density and writes `optimized_flow.png` and
`optimized_flow.npz` in the same results folder.

Pressure is fixed to 0 Pa gauge at both outlets; the flow split is a solver
result, not a prescribed ratio. The first mesh is coarse and the discretization
uses an orthogonal-grid gradient approximation. Treat results as a solver
prototype until residuals, global mass balance, and mesh sensitivity have all
been checked.

Each log line reports local continuity, overall boundary mass balance, field
change, and the linear-system residual. A run is marked converged only when all
four are below the configured tolerance. The outlet split and coarse-mesh
solution still need comparison with a trusted benchmark before design use.

## Field definitions

- `u`, `v`: mean x and y velocity (m/s); `p`: gauge pressure (Pa).
- `k`: turbulent kinetic energy (m^2/s^2); `omega`: specific dissipation
  (1/s); `nut`: turbulent kinematic viscosity (m^2/s).
- `F1` and `F2` are the SST wall/blending functions. `nut` uses the SST shear
  stress limiter. The k and omega equations use blended constants, production
  limiting, and omega cross diffusion.
- `mass` is the normalized discrete continuity residual after face-flux
  correction; `change` is the largest normalized field update. Both should
  decrease before interpreting the saved fields.
- `pseudo_cfl` controls the stabilizing pseudo-time term; `velocity_relaxation`,
  `pressure_relaxation`, and `velocity_correction_relaxation` control iteration
  damping in `case.py`.
