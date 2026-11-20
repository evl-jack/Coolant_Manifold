# 2D Brinkman–RANS verification and optimization

## Models

The density convention is `gamma=1` for fluid and `gamma=0` for solid. Model A adds
the Darcy/Brinkman momentum sink `-alpha*u`, where
`alpha=alpha_max*(1-gamma)**power`. Models B and C use the same penalty to suppress
solid flow and reconstruct the `gamma=0.5` interface with linear or guarded
quadratic edge interpolation. They use the reconstructed wall distance in the
SST closure and a near-interface momentum damping term as an immersed no-slip
approximation. Model A retains the baseline rectangular-wall distance.

Run the original rectangle (`python -m rans_solver.solve --case rectangle`) and
the explicit cut-cell T reference (`--case manifold`) with `--model explicit`.
Run density models with `--model A`, `--model B`, or `--model C`. For a four-way
CSV comparison, run `python -m rans_solver.compare_models`.

## Optimization

`python -m rans_solver.optimize --steps 5 --model B` runs a deliberately coarse
finite-difference demonstration. It uses a 4-by-3 design around both main
channel walls in the middle 15 mm of the domain. The first and last 2.5 mm
of the main channel are fixed. A fixed 2.25 mm wide branch root spans
y = 6.00–6.25 mm and steps into the 1.25 mm wide branch. The
inlet and main outlet span y = 2.25–5.75 mm. Each opening ends in a split
boundary face with a no-slip wall segment. The density mesh starts at y = 0 mm,
giving the lower wall 2 mm of design space. From y = 8.0 to 9.5 mm, the side
branch is fixed at 1.25 mm wide and the surrounding cells are fixed as solid.
The run starts from the T density. It keeps a central main-channel spine and
a branch-to-outlet spine fixed, checks
connectivity, and fixes the fluid-equivalent area at 30% of the editable region
by default. Set `--volume-fraction` to override that constraint. The side branch
is 1.25 mm wide. The optimizer targets `Q_side = 0.75 * Q_main`, or about
42.86% of total discharge through the side outlet. The constraint allows an
absolute mismatch of 0.1% of inlet mass flow in that equation and requires outward
flow at both outlets. Set `--side-to-main-ratio` to change the target.
Both outlets use pressure boundary conditions. Candidate
solves warm-start from the base solution and must meet a looser `5e-4` field
change tolerance using up to 1000 flow iterations. The search reduces the
outlet balance error until it finds a feasible design, then accepts only feasible
steps that lower pressure drop. If no feasible design is found, it saves the
lowest-violation iterate and records its status in `optimization.csv`.
Brinkman resistance defaults to `1e6` in the optimizer. Results are written
under `rans_solver/results/`.

### Search using channel dimensions

Run `python -m rans_solver.shape_search --steps 5` to vary three lengths in
millimeters: the main passage narrowing before and after the junction, and
the flare of the side branch below its fixed outlet section. The smooth wall
profiles keep one opening per cross section. They preserve the fixed sections,
keep the main passage at least 2 mm wide, and use the 30% editable volume
setting as a maximum. Cell densities for the search come from geometry
coverage and stay between zero and one.

The search first tries to match the 3:4 outlet flow ratio. After reaching the
ratio tolerance, it can reduce the static pressure drop while keeping the
ratio and volume constraints. The gradient CSV compares finite differences
at 0.08 and 0.04 mm. The candidate CSV records convergence, geometry and
volume checks, and flow results. The explicit wall CSV checks the binary
geometry at 0.5 and 0.25 mm for convergence, outlet balance, volume, and wall
leakage. A nearby gradient or a locally best candidate does not establish
whether every possible pipe shape can meet the target.

## Interpretation and limits

The explicit case is the reference geometry; compare its mass split, pressure
drop, and conservation against A/B/C at matching resolution and boundary data.
The solver is an educational first-order finite-volume SIMPLE implementation,
not a validated engineering package. The current immersed condition is a
near-wall damping approximation rather than a formally conservative sharp
immersed-boundary flux treatment. The optimizer uses finite differences rather
than an adjoint, and its objective/volume projection are starter choices. Check
mesh refinement, penalty sensitivity, wall reconstruction, and sensitivity
convergence before using any design result. No 3D case is included.
