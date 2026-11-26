# Coolant Manifold Project 4 Handoff

**Status on November 20, 2026.** The self-written 2D finite-element solver, its benchmark adjoint, and a fixed-geometry T-junction baseline are working on ASU Sol. The next development milestone is a verified fillet-only design formulation. The existing finite-volume RANS and Brinkman studies remain the research record.

The subsequent fillet study verified the 20-control design and manifold gradients. The tested cases did not reach the proposed flow band. See the current [formulation and feasibility decision](FORMULATION.md).

The course deliverable calls for verified in-house FEM, checked adjoint gradients, constrained optimization with KKT analysis, a comparison of optimized designs with and without nonlinear physics, and a GitHub Pages report with a live initial-versus-optimized WebGL view. The course page lists November 25 for presentations and the final site. A 2D passage extruded for display meets the 3D geometry requirement. The pasted roadmap ends partway through the design-variable section, so the final control layout and engineering limits are still open. [Official Project 4 specification](https://designinformaticslab.github.io/DesignOptimization2025/project4.html)

## Code and saved evidence

| Item | Role |
|---|---|
| [FEM solver](fem.py) | Taylor–Hood P2 velocity and P1 pressure on affine triangles. Includes viscous stress, pressure and continuity coupling, Brinkman drag, loads, traction, nonlinear convection, a damped Newton solve, and a drag adjoint. |
| [Verification driver](verify.py) | Runs the nine benchmark checks and writes measured errors, thresholds, and pass flags. |
| [Fixed T-junction driver](junction.py) | Builds two boundary-fitted meshes from the saved geometry dimensions and runs Stokes, Re 25, and Re 50 by continuation. |
| [Latest verification report](output/sol-2082615/verification.json) | All nine checks passed in Sol job **2082615**. |
| [Fixed T-junction report](output/sol-2082569/junction-reference.json) | Numerical results from Sol job **2082569**. The [coarse state](output/sol-2082569/reference-grid-1.npz) and [fine state](output/sol-2082569/reference-grid-2.npz) contain mesh coordinates, connectivity, and the final Re 50 state. |
| [Geometry checkpoint](../output/max-side-flow-sol-2079305/summary.json) | Original Sol sweep record. The FEM driver reads the four geometry dimensions from its `final` entry. |

The state vector stores all horizontal velocities, all vertical velocities, then pressures. The mixed formulation uses symmetric-gradient viscosity. All-velocity benchmarks pin one pressure degree of freedom. In the junction case, natural zero traction at both outlets establishes the pressure datum. The adjoint solves the transpose of the converged free-degree Newton tangent and contracts the adjoint and forward velocities with the resistance derivative. The current adjoint covers fixed mesh and boundary data with design dependence entering through resistance. A response that depends directly on design also needs its explicit derivative.

## Verification and baseline

| Check | Observed result |
|---|---|
| Basis and quadrature | Nodal interpolation, partition of unity, gradient-sum, and monomial integration checks passed. |
| Analytical Poiseuille flow | Relative velocity \(L_2\) error \(1.18\times10^{-14}\). |
| Invariance and conservation | The rotated channel reproduced the exact solution. Tested channel mass balances passed. |
| Manufactured Stokes refinement | Observed velocity \(L_2\) order **3.10** between the finest two grids. |
| Manufactured nonlinear flow | The residual and Newton tangent check passed. Tangent relative error was **\(1.08\times10^{-11}\)**. |
| Brinkman drag and leakage | The constant-drag benchmark passed. In the layered benchmark, the fraction of outlet flow through nominal solids fell from **21.12%** to **4.11%** as nondimensional resistance rose from 100 to 1,600. |
| Drag adjoint | Two channel responses and 24 controls. At the smallest tested step, the largest directional relative error was below **\(1.4\times10^{-7}\)**. Taylor remainder orders were close to two. |

The junction uses a **1.000 mm** main-bore radius, **0.625 mm** side-outlet radius, **0.500 mm** branch-throat radius, and **1.500 mm** junction fillet radius. The course case uses a parabolic inlet, unit nondimensional inlet flow per depth, no-slip walls, and equal zero-traction outlets. The length and speed scales are the 2 mm main-bore height and mean inlet velocity. **Re 50** is the selected reduced laminar operating point.

| Mesh | Stokes side share | Re 50 side share |
|---|---:|---:|
| 708 triangles | 27.389% | 24.990% |
| 2,832 triangles | 27.399% | 24.905% |

On the fine grid, the Re 50 side/main ratio was **0.33165** and the relative mass imbalance was **\(2.66\times10^{-15}\)**. The two Re 50 meshes differed by **0.085 percentage points** of inlet flow. This is a two-grid sensitivity measurement for the fixed geometry. The program has benchmark adjoint evidence and still needs gradient checks for the eventual manifold responses.

## Relationship to the earlier RANS work

The research checkpoint reached **0.1852 side-to-main** and **15.63% side share** with explicit walls. The fresh coarse re-solve reproduced that result. Two fine-grid RANS attempts ended near **18%** and missed the steady momentum and pressure-correction gates. The earlier reconstructed contour selected the same whole cells as the reference on each grid, and that study stopped before measuring Brinkman hydraulic differences. The laminar FEM baseline uses its own inlet profile, operating point, and outlet formulation. An independent CFD comparison should match the governing physics and boundary conditions first. [RANS verification](../output/junction-verification-sol-2081608/reference/reference.json) · [Reconstruction audit](../output/junction-verification-sol-2081608/retry/geometry-audit.json)

The research budget counted **2.211 mm²** above the old 6 mm datum, while the complete editable branch from the current bore wall at 5 mm measured **4.164 mm²** on that coarse mesh. Define the FEM budget against the full intended editable region, including the fillet, before comparing candidate areas.

## Reproduction

From the workspace root in a Sol compute allocation, use the existing scientific Python environment with NumPy and SciPy. The commands write fresh results to `project4_fem/output/`. The archived results are under job-specific subfolders.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -m project4_fem.verify
python -m project4_fem.junction \
  --checkpoint output/max-side-flow-sol-2079305/summary.json \
  --reynolds 50
```

The verification command exits nonzero if any declared check fails. The original Sol source and outputs were staged under `/scratch/$USER/manifold-rans-20261119/project4-fem.s6BLfS`.

## Next milestone

1. **Freeze the course design problem.** Define at least 20 independent controls in one fillet region and mirror their effect to the other. Keep port widths, locations, and straight runs fixed. Choose the side-flow requirement, a distinct second PDE-dependent constraint, and the complete fillet area limit. Check feasibility at Re 50. If a 0.25 side/main target is retained, it corresponds to **20% side share**, below the current laminar baseline.
2. **Verify the fillet design map.** Use a smooth density field on a fixed FEM mesh. Measure leakage near the fillet as the penalty and interface resolution vary. Check the reconstructed geometry for symmetry, connectivity, and minimum width. Run the matched independent laminar CFD comparison when its boundary conditions are defined.
3. **Verify the manifold gradients.** Compare the adjoint of each final response with finite differences and Taylor tests on converged manifold states. Keep thresholds, perturbation sizes, and elapsed time in the results.
4. **Optimize and report.** Run the same design problem in Stokes and Navier–Stokes, verify the selected geometry in a fresh explicit-wall solve, and report feasibility, KKT stationarity, complementarity, active constraints, and multipliers. Measure optimization runtime and generate the GitHub Pages report and WebGL extrusion from the saved initial and optimized 2D designs.
