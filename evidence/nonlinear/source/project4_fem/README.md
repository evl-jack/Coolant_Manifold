# Project 4 finite-element model

This module implements the reduced 2D flow analysis for the coolant manifold. It shares the frozen T-junction dimensions with the earlier Sol study. The course model uses laminar flow at its own nondimensional operating point. The prior RANS results remain evidence for the finite-volume research model.

## Analysis

`fem.py` assembles continuous quadratic velocity and linear pressure on affine triangles. It implements the symmetric-gradient viscous term, mixed pressure and continuity blocks, Brinkman drag, body forces, and natural traction. The nonlinear residual includes Navier–Stokes convection. A sparse direct solve supplies each damped Newton step. The Newton tangent includes derivatives with respect to both velocity factors in convection.

The weak continuity equation is multiplied by −1 to give a symmetric Stokes saddle-point matrix. All-velocity benchmarks pin one pressure degree of freedom. The pressure datum in the T-junction is set by the two zero-traction outlet boundaries. Fluxes use exact Simpson integration of the quadratic normal velocity along the boundary edges.

For a converged state, `drag_adjoint` solves the transpose of the free-degree tangent for each response. Each gradient contracts the forward velocity and adjoint velocity with the quadrature-point derivative of the drag. This interface covers fixed boundary data and design variables that enter only through resistance. A directly design-dependent response also needs its explicit design derivative.

## Reproduce the Sol checks

Use the existing Sol scientific Python environment from the project root. Set `OPENBLAS_NUM_THREADS=1` and run these commands in a Slurm allocation.

```bash
python -m project4_fem.verify
python -m project4_fem.junction --checkpoint output/max-side-flow-sol-2079305/summary.json
```

The verification script reports the element and quadrature checks, analytical channel flow, coordinate rotation, global mass conservation, manufactured-solution mesh convergence, a manufactured nonlinear flow with a tangent check, analytical Brinkman channel flow, leakage through Brinkman layers, and directional adjoint checks with Taylor orders. It saves the measured values and thresholds in `project4_fem/output/verification.json` and returns a nonzero exit status if a check fails.

The T-junction driver solves the fixed geometry on two boundary-fitted meshes. Geometry edges are piecewise linear and every triangle uses the same mixed element. The inlet has a parabolic velocity profile of unit mean velocity. The nondimensional scales are the main channel height and inlet mean velocity. The default Reynolds number is 50 and the driver uses Stokes and half-Reynolds solutions to initialize the final Newton solve. It saves both resolutions and their state and flow data. The reported two-grid difference measures sensitivity to those grids. A matched independent CFD comparison is a separate verification rung.

## Sol results

All nine entries in the verification ladder passed. The manufactured velocity error decreased at order 3.10 on the finest grid pair. Two response gradients checked with 24 resistance controls had directional relative errors below `1.4e-7` at the smallest tested step and Taylor remainder orders near two. The layered channel tests measured less solid-region leakage as the Brinkman resistance increased. [Verification results](output/sol-2082615/verification.json)

At Reynolds 50 the frozen junction sent 24.99% of inlet flow to the side on 708 triangles and 24.91% on 2,832 triangles. The finer mesh gave 27.40% in Stokes flow. These values establish a baseline for the reduced course model at the stated boundary conditions. [Junction reference](output/sol-2082569/junction-reference.json)

## Course design phase

1. Define at least 20 independent mirrored design parameters inside the two fillet regions. Keep the straight runs and all port dimensions fixed. Use a smooth density map on a fixed FEM mesh so the discrete adjoint can differentiate the resistance directly.
2. Declare the optimization objective and at least two flow-dependent constraints with explicit limits. Check feasibility near the baseline before committing to a side-to-main target. Measure the complete editable branch and fillet area from the current bore wall when defining the area cap.
3. Check every final response gradient with finite differences and a Taylor test on converged manifold states. Then compare a Stokes optimum with a Navier–Stokes optimum at the same operating point and report the KKT residuals, active limits, multipliers, and runtime.
4. Build the report and interactive initial-versus-optimized extrusion from the saved 2D results. Keep the 2D performance analysis and the 3D display geometry traceable to the same design files.

[Official Project 4 specification](https://designinformaticslab.github.io/DesignOptimization2025/project4.html)
