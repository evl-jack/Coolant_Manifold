# Coolant manifold RANS project handoff

Status through November 19, 2026. This is a chronological record of the mesh, solver, optimization, and ASU Sol work. Values from intermediate runs are approximate where indicated. Current CSV files cover the latest recorded runs and may supersede older outputs.

## 1. Cartesian mesh and first 2D verification

- Started in `mesh.py` with a small-grid connectivity check for fluid cells, then simplified the mesh code and inspected the grid visually. The mesh is Cartesian. The manifold mesh adds grid lines at geometric landmarks, so cells can have unequal spacings. Early setup and plots are in `analysis.py`, `ANALYSIS.md`, and `output/`.
- Reviewed Lars Davidson’s pyCALC-RANS PDF and extracted the supplied `pythons-rans-code-RANS-open.tar.gz` into `pyCALC-RANS-reference/`. The plan changed from using that package directly to writing a concise in-house 2D solver in `rans_solver/`. The package was retained as a methodology reference. [Davidson’s pyCALC-RANS notes](https://www.cfd-sweden.se/lada/postscript_files/py-calc-rans.pdf)
- Made the first run a 20 × 4 mm rectangle with no side outlet. The nominal 1 mm grid gives 80 cells. Water properties were $\rho=998.2$ kg/m³ and $\mu=1.002\times10^{-3}$ Pa·s. The inlet is 2 m/s along x, the outlet is 0 Pa gauge, and the inlet turbulence inputs are assumed 5% intensity and 0.56 mm length scale.
- `rans_solver/solve.py` uses cell-centered finite volumes, SIMPLE outer iterations, Rhie–Chow face fluxes, first-order upwind convection, and k–ω SST. SciPy `spsolve` directly solves each sparse transport and pressure-correction system. Pseudo-time stabilization and under-relaxation improve outer convergence. The rectangle converged in about 102 iterations. The first explicit-wall T added a 1 mm branch centered at x = 10 mm from y = 4 to 7 mm, and converged in about 420 iterations on the coarse mesh.
- Reviewed pressure, velocity magnitude, and velocity direction plots. Direction comes from the mean velocity components $(u,v)$, shown with arrows. The pressure gradient belongs in the momentum equation. Convergence requires local continuity, global mass balance, field change, and linear residual to pass the applicable tolerance.

## 2. Brinkman coupling and first density optimization

- Implemented the supplied Brinkman–RANS plan in 2D. The 3D verification was deferred at this stage. $\gamma=1$ is fluid and $\gamma=0$ is solid. Model A adds the implicit momentum sink $-\alpha \mathbf{u}$ with $\alpha=\alpha_{\max}(1-\gamma)^2$. Models B and C reconstruct the $\gamma=0.5$ contour with linear or guarded quadratic interpolation and use the reconstructed wall distance for SST and near-wall damping. `brinkman.py`, `compare_models.py`, and `BRINKMAN.md` contain the implementation and comparison setup.
- The $10^8$ maximum penalty was stiff in early density runs. The optimization path uses $10^6$, warm starts and a $5\times10^{-4}$ tolerance with up to 1,000 iterations per solve. The baseline default is $10^{-5}$. The later shape search uses $10^{-4}$. The B/C immersed wall treatment is approximate. Explicit wall boundary faces enforce no-slip and zero normal flux.
- Improved `optimize.py` after early runs barely changed the geometry. The revised search uses a 4 × 3 coarse control field over the fine density grid, area-weighted volume projection, fixed ports and channel spines, four-neighbor connectivity at the chosen density threshold, feasible finite-difference perturbations, convergence checks, and backtracking. Saved density is rerun by `plot_optimized.py`, which plots contours and velocity vectors on the nonuniform grid.
- Early objective used pressure loss plus a soft 1% inlet-flow side target. A short 4 × 3 run reduced the objective from about 388 to 378, with side share still around 4.15%. The first upper-band edit mask had 34 cells. An 8% editable-fluid budget mostly added sub-threshold density. At 20%, 18 of the 34 cells moved by more than 0.05 in density, around three crossed the $\gamma=0.5$ boundary, and the side share was around 4.7%.

## 3. Design region, outlet targets, and volume revisions

In the table, “side share” is $Q_\mathrm{side}/Q_\mathrm{in}$. “Ratio” is $Q_\mathrm{side}/Q_\mathrm{main}$. The reported intermediate values come from different geometries and objectives.

| Order | Change and execution | Observed outcome |
|---|---|---|
| 1 | Widened the editable region around fixed inlet and main-outlet planes at x = 0 and 20 mm. Raised the side outlet by 0.5 mm from y = 7.0 to 7.5 mm. | 97 editable cells in one run, about 7.5% side share. The lower main wall still sat on the outer boundary at y = 0. |
| 2 | Moved the whole T up 2 mm inside the mesh. Main passage moved to y = 2–6 mm and the side outlet to y = 9.5 mm. The enclosure floor stayed y = 0, giving the bottom wall 2 mm of design space. | 22 × 10 density grid, 136 editable cells at that stage. The threshold shape barely changed. Side share went from about 5.41% to 5.14%. |
| 3 | Widened the side opening by 0.25 mm total, from 1.0 to 1.25 mm, centered at x = 10 mm. Introduced a hard 1% side-share cap and ranked candidates by violation before pressure drop. | Side share fell from about 6.19% to 4.95%. Pressure drop went from about 319 to 777 Pa. The cap was unmet. |
| 4 | Replaced the 1% cap with $Q_\mathrm{side}=0.75Q_\mathrm{main}$. Both exits remained at 0 Pa gauge. | Target corresponds to 42.86% side share for balanced flow. After five steps with the original 4 mm main openings, side share reached about 23.03%, ratio 0.299. |
| 5 | Trimmed each main opening by 0.25 mm on both edges about y = 4 mm. New opening spans y = 2.25–5.75 mm, height 3.5 mm. Split boundary faces preserved the existing density grid and added wall segments outside each opening. | Inlet flow per unit depth dropped 12.5% to 0.007 m³/s per meter. Side share was about 23.12%, ratio 0.301. Sampled speed along the B contour still reached about 1.05 m/s, down from about 1.21 m/s before the trim. |
| 6 | Froze 2.5 mm at both ends of the main passage. Widened the branch root by 0.5 mm on each side for its first 0.25 mm in y. | Current root is 2.25 mm wide from x = 8.875–11.125 mm and y = 6–6.25 mm. The remaining side bore is 1.25 mm. Side share was about 23.16%, ratio 0.301. |
| 7 | Investigated a bulge below the side outlet. A threshold cell around $\gamma=0.516$ linked the spur near y = 8.5 mm to the passage. Froze y = 8–9.5 mm to the reference bore and surrounding solid. | The opening stayed 1.25 mm. The bulge moved lower to around y = 7.5 mm. Side share decreased to about 22.17%, ratio 0.285. |
| 8 | Reduced the editable-area fluid budget from the reference equivalent of 34.73% to 30.0%. | The lower spur shrank to three threshold cells covering 2.75 mm². The latest density run reached 22.61% side share, ratio 0.2922, with about 653 Pa pressure drop. |

Current 2D constraints live in `mesh.py` and `rans_solver/case.py`. The editable main span is between x = 2.5 and 17.5 mm, with paths and ports protected. The reference fluid area is 84.625 mm². Of 87.71875 mm² editable area, the reference uses 30.46875 mm² fluid-equivalent. The new 30% budget gives 26.315625 mm² editable fluid-equivalent, or 80.471875 mm² total including fixed fluid. In the density search the budget is enforced as an equality using $\sum_i A_i\gamma_i$, with $0\leq\gamma_i\leq1$.

With the trimmed inlet, the target split is 0.004 main and 0.003 side m³/s per meter of depth. The optimizer requires outward flow at both outlets and uses $\left|Q_\mathrm{side}-0.75Q_\mathrm{main}\right|/Q_\mathrm{in}\leq0.001$. Outlet pressure boundary conditions and the geometry determine the actual split.

## 4. Gradient-based passage search and 2D wall checks

The next plan moved from editing coarse density cells to optimizing continuous passage profiles. `rans_solver/shape_search.py` uses three geometric controls in mm, each measured per wall or edge. The upstream and downstream main-passage pinches range from 0 to 1.0 mm. A side-passage flare ranges from 0 to 0.5 mm. Smooth sine-squared profiles preserve the ports and fixed segments and keep each passage a single open interval. Minimum main width is 2 mm. Cell coverage supplies Brinkman density, and both sampled and analytic areas must stay within the 80.471875 mm² cap.

SLSQP uses finite-difference gradients to reduce outlet imbalance, with pressure loss secondary once the flow constraint is feasible. Baseline side share was about 12.71% for this search. Finite differences agreed on the local trends at two perturbation sizes. The downstream pinch had the strongest sensitivity.

| Control | Side-share sensitivity at 0.08 mm, percentage points/mm | Side-share sensitivity at 0.04 mm, percentage points/mm |
|---|---:|---:|
| Upstream pinch | −0.780 | −0.814 |
| Downstream pinch | +14.424 | +14.540 |
| Side flare | +0.344 | +0.340 |

The five-step search evaluated 44 candidates including probes. Best controls were **0.003734 mm upstream pinch, 1.0 mm downstream pinch, and 0.268893 mm side flare**. The selected shape uses 78.6968 mm² analytic area and 78.6888 mm² on the density grid. Tested it with explicit wall boundary faces on binary masks formed from 50% cell coverage.

| Model or spacing | Side share | Side/main | Pressure drop |
|---|---:|---:|---:|
| Coarse Brinkman B search model | 20.57% | 0.2590 | 594 Pa |
| Explicit walls, 0.50 mm | 28.06% | 0.3900 | 1,784 Pa |
| Explicit walls, 0.25 mm | 28.60% | 0.4005 | 1,520 Pa |

Both explicit runs converged, with zero modeled wall flux and inlet/outlet mass balance to numerical precision. At 0.25 mm the explicit flow plot took 1,327 iterations. The code produces a working 2D solution, but even the finer explicit result misses the 0.75 target. The search was limited to these three controls, and pressure remains mesh-sensitive. There is no independent validation or proof that a feasible shape is impossible.

Use `rans_solver/results/shape_flow.png` for the selected shape’s explicit-wall flow plot. `optimized_flow.png` is the earlier density design plot. `shape_best.csv`, `shape_candidates.csv`, `shape_gradients.csv`, and `shape_wall_checks.csv` record the search. `optimization.csv` records the density history for the latest density run.

## 5. 3D geometry, mesh, and Sol execution

- Reopened 3D verification after the 2D work. The user first requested a plan only, then directed the case to be built in VS Code and run on ASU Sol. `sol_rans_3d/build_mesh.py` reads `shape_best.csv`, revolves the main passage around x through y = 4 mm and z = 0, revolves the branch around y through x = 10 mm and z = 0, and unions both volumes into one fluid region. This gives a single circular side outlet. Main openings are 3.5 mm diameter with 0.25 mm axial lips. The downstream throat is 2 mm, the side exit is 1.25 mm, and the branch root is 2.25 mm. CAD fluid volume is about 223.54 mm³. The 2D area constraint did not set a 3D volume limit.
- Gmsh generated tagged meshes in meters and a STEP file in millimeters. Main opening area is 9.621 mm² each and side opening area is 1.227 mm². Coarse mesh has 219,147 tetrahedra with nominal 0.27 mm core and 0.075 mm wall sizing. Fine mesh has 622,031 with 0.19 mm core and 0.05 mm wall sizing. Gmsh output showed positive shape metrics, with minimum SICN about 0.0345 and 0.0184 respectively. OpenFOAM performs its own checks.
- `sol_rans_3d/case/` configures steady incompressible OpenFOAM `simpleFoam`, SST with wall functions, water, 2 m/s inlet, zero gauge pressure at both outlets, and no-slip walls. Inlet $k=0.015$ m²/s² and $\omega=399.30$ s⁻¹. The prescribed inlet corresponds to about 1.155 L/min. `submit.sbatch` requests one node, 16 MPI ranks, 32 GB and four hours on `public`. `run_openfoam.sh` imports the tagged mesh, assigns the wall patch type, and gates the solve on `checkMesh -allGeometry -allTopology`. `summarize.py` reports outlet flux, mass balance, pressure drop and other checks after the solve. The job writes `manifold.foam` for ParaView.
- The plan initially used an OpenCFD v2312 Apptainer image. After SSH access worked, the actual jobs used Sol’s native OpenFOAM v2512, with a shell initialization fallback in `submit.sbatch`. The browser shell was available for the user’s login and key registration. Host identity was checked against the authenticated Sol shell. After clearing a public-key clipboard issue, SSH access and the uploaded bundle checksum were verified. The remote case lives at `/scratch/$USER/manifold-rans-20261119/sol_rans_3d` on Sol. The runner with the dual-mesh option was uploaded separately after unpacking the original bundle.
- **Job 2053128 failed the coarse mesh gate before the solver started.** `checkMesh` reported 700 cells with determinant below 0.001 and 14 faces with interpolation weight below 0.05. The minimum weight was about 0.0182. Maximum nonorthogonality was about 81.1°, with four faces above 70°. The failure occurred after about 15 seconds. No flow fields or flow split resulted.
- **Job 2057077 failed the dual-mesh retry.** It used `polyDualMesh -overwrite -concaveMultiCells 45`. The low-determinant count fell from 700 to 22, but the mesh failed eight checks, including face orientation and severe skewness. The solver did not start.

## 6. Hex-mesh repair and Sol restart

- In the original tetrahedral mesh, 682 of the 700 low-determinant cells had only two internal faces. Of those 700, 534 were around the main ends or lips and 166 around the branch. All 14 low-interpolation-weight faces were around the branch.
- Exported the same CAD skin as a closed, four-region STL and generated hex-dominant grids with `snappyHexMesh`. Local refinement removed the original failures but left concave cells at internal refinement transitions. Uniform grids reduced the residual issues to 23 cells on the 0.09 mm coarse grid and 12 on the 0.06 mm fine grid. Decomposing cells and collapsing sliver edges introduced additional problems.
- `repair_hex_mesh.py` uses a bounded nonlinear least-squares fit to adjust vertices around the flagged cells. It minimizes violations of face-plane convexity and keeps the port vertices fixed. Maximum wall point motion was 0.01455 mm on the coarse grid and 0.00957 mm on the fine grid. Both passed the unchanged `checkMesh -allGeometry -allTopology` gate on Sol. Coarse job 2059820 verified 306,218 cells and minimum interpolation weight 0.1215. Fine job 2059943 verified 1,036,985 cells and minimum interpolation weight 0.1060. Both have one connected region. The meshed side outlet area is 1.18% below CAD on the coarse grid and 0.49% below on the fine grid.
- The separate Sol working directory is `/scratch/$USER/manifold-rans-20261119/uniform-mesh`. The repaired `polyMesh` directories are inside `repaired-coarse/` and `repaired-fine/`. `MESH_PREBUILT` lets the job reuse these grids and runs the same geometry and topology gate. Coarse solver job 2059887 passed the mesh gate and decomposition, then hit an Open MPI CPU-binding error. The runner was switched to Slurm's `srun --mpi=pmix` for native OpenFOAM. Coarse solver retry 2059974 was running at iteration 376 at the last check. [ASU MPI launcher instructions](https://docs.rc.asu.edu/mpi/)

## 7. Where to resume

Use `sol_rans_3d/README.md` for setup. From an authorized Sol session, check `squeue -j 2059974` and inspect `runs/coarse-2059974/log.simpleFoam`. Confirm solver convergence, then read `summary.json` and inspect the reconstructed fields. The fine mesh is ready for a separate solve. Compare pressure drop, outlet fluxes, global balance, and wall $y^+$ between the two resolutions. The main design question remains whether the side/main ratio can reach 0.75 at the current pressure boundary conditions, volume budget, and geometry constraints.

