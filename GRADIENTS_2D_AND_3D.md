# How we used gradients in the coolant manifold

In **2D**, we measured how changing the passage shape on a Cartesian grid affected flow split and pressure drop. In **3D**, we measured how moving mesh vertices affected the quality of nearby cells. The RANS solver also computes spatial gradients such as $\nabla p$ inside each flow solve.

## 2D · Optimize the passage on a Cartesian grid

### Build and check the grid

The first mesh fix was geometric. We placed grid lines at the ports, branch root, and fixed channel sections and omitted nearby regular lines that would create sliver cells. On a small grid, a graph traversal checked that every fluid cell was reachable through shared faces. The test also checked the boundary tags, port lengths, and reference fluid area. [Grid generation and small-grid check](mesh.py)

### Use gradients to change which cells carry fluid

The first optimizer perturbed a **4 × 3** array of controls and interpolated the changes onto the density grid. A density $\gamma$ ranges from zero for solid to one for fluid. Brinkman Model B suppresses flow in solid and gray cells. Before each trial solve, the optimizer projected the densities back to the area-weighted fluid budget, preserved fixed regions, and checked that the inlet connected to both outlets. It used finite differences and backtracking to accept steps that reduced the outlet mismatch, with pressure drop as a secondary criterion. This could change a cell’s fluid fraction by a small amount without changing the visible $\gamma=0.5$ boundary. [Density optimizer](rans_solver/optimize.py) · [Brinkman models](rans_solver/BRINKMAN.md)

We then switched to three **geometric controls** in millimeters. They define smooth sine-squared changes to the width while preserving the inlet, both outlets, and the fixed branch root. The main pinches range from 0 to 1.0 mm per wall and the side flare from 0 to 0.5 mm per edge. Each passage occupies a single open interval across its width. Cell coverage maps the shape back onto the Cartesian grid. The area cap was **80.471875 mm²** including fixed fluid, with at most 30% of the editable area allocated to fluid. [Shape and area constraints](rans_solver/shape_search.py)

For each control $a_i$, we perturbed it in either direction, ran RANS for both candidates, and measured the change in a flow metric $F$. When a perturbation hit a control bound, we used the actual distance between the sample points.

$$
\frac{\partial F}{\partial a_i}
\approx
\frac{F(a_i^+)-F(a_i^-)}{a_i^+-a_i^-}
$$

Both outlets were set to zero gauge pressure. The design target was $Q_\mathrm{side}/Q_\mathrm{main}=0.75$, equivalent to a 42.86% side share for balanced flow. The bounded SLSQP search minimized $\tfrac12 e^2$, where $e=(Q_\mathrm{side}-0.75Q_\mathrm{main})/Q_\mathrm{in}$. Its gradient is $e\nabla e$. The pressure-drop stage was gated on first meeting the flow mismatch tolerance of 0.001 relative to inlet flow. The tested candidates never reached that tolerance. [Flow metrics and tolerance check](rans_solver/optimize.py) · [Search implementation](rans_solver/shape_search.py)

| Control | Side-share change at 0.08 mm perturbation | Side-share change at 0.04 mm perturbation |
|---|---:|---:|
| Upstream pinch | −0.78 percentage points/mm | −0.81 percentage points/mm |
| Downstream pinch | +14.42 percentage points/mm | +14.54 percentage points/mm |
| Side flare | +0.34 percentage points/mm | +0.34 percentage points/mm |

The downstream pinch had the largest local effect at the starting design. The five-step search evaluated 44 candidates including gradient probes. The selected shape used a **0.0037 mm upstream pinch**, **1.0 mm downstream pinch**, and **0.2689 mm branch flare**. Its analytic area was **78.70 mm²**. We rebuilt the shape as a binary fluid mask at 50% cell coverage and solved it with explicit no-slip faces. The 0.50 and 0.25 mm runs both converged with zero modeled wall flux. Their side/main ratios were **0.3900** and **0.4005**, below the 0.75 target. The geometry passed the connectivity check and became the input to the 3D CAD build. [Measured gradients](rans_solver/results/shape_gradients.csv) · [Selected design](rans_solver/results/shape_best.csv) · [Explicit-wall checks](rans_solver/results/shape_wall_checks.csv)

## 3D · Repair the mesh around the flagged cells

We built the 3D fluid volume by revolving and joining the selected 2D passage profiles. The original tetrahedral mesh failed OpenFOAM’s cell determinant and face interpolation checks. We rebuilt the mesh with uniform hex-dominant grids. The replacement passed those two checks and left **23 coarse** and **12 fine** cells flagged for concavity. These remaining defects were the starting point for the gradient repair. [3D mesh setup](sol_rans_3d/README.md)

For each flagged cell and the neighbors affected by moving its vertices, the repair evaluates where the other face centers sit relative to each outward-facing plane. If $\mathbf c_f$ is the center of face $f$ and $\mathbf n_f$ its outward unit normal, then a positive value of the following expression means face center $\mathbf c_g$ is on the exterior side of that plane.

$$
d_{fg}=(\mathbf c_g-\mathbf c_f)\cdot\mathbf n_f
$$

The fit penalizes face centers that cross the plane or come too close to it. It also penalizes faces that sit too close to an estimated cell center and discourages unnecessary point motion. The variables are the $x,y,z$ displacements of nearby mesh vertices, scaled by nominal cell spacing. Each residual depends on the vertices of one cell, so SciPy estimates a **sparse finite-difference Jacobian** and uses bounded nonlinear least squares to reduce the combined residuals. The inlet and outlet vertices stay fixed. Wall vertices have tighter motion bounds than interior vertices. The repair writes a replacement `points` file and preserves mesh connectivity. [Vertex repair implementation](sol_rans_3d/repair_hex_mesh.py)

| Grid | Final cell count | Flagged concave cells before → after | Max wall point motion | Side outlet area below CAD |
|---|---:|---:|---:|---:|
| 0.09 mm coarse | 306,218 | 23 → 0 | 0.01455 mm | 1.18% |
| 0.06 mm fine | 1,036,985 | 12 → 0 | 0.00957 mm | 0.49% |

Both repaired grids passed `checkMesh -allGeometry -allTopology` on Sol as single connected fluid regions. The RANS runner reused the repaired coarse mesh and ran the same quality gate before starting the solver. Mesh acceptance and flow convergence are separate checks. [Validation jobs and solver handoff](MANIFOLD_HANDOFF.md#6-hex-mesh-repair-and-sol-restart) · [OpenFOAM runner](sol_rans_3d/run_openfoam.sh)

## Reproduce the gradient steps

Run the 2D search from the workspace root. It rewrites the shape result CSVs.

```powershell
.\.venv\Scripts\python.exe -m rans_solver.shape_search --steps 5
```

For 3D, copy the original `polyMesh` directory with its `sets/concaveCells` file. Use `0.09` for coarse or `0.06` for fine. Put the output in a copy of the mesh as `points` and run the full OpenFOAM mesh check on that copy.

```powershell
.\.venv\Scripts\python.exe sol_rans_3d/repair_hex_mesh.py .\polyMesh 0.09 .\points.repaired
```

