# SST closure wall-distance revision

The frozen binary T uses Euclidean distance from each quadrature point to the nearest of the five physical wall segments. The open port planes are excluded. The revised run uses `d_eff = max(d_geom, c h_local)` for SST blending, with `h_local = 2 A / longest_edge`, the minimum triangle altitude. Both lengths use inlet height as the reference. The Sol run sets `c=0.5`; the driver default of zero retains the original distance for replay. This floor is a numerical modeling choice.

The scalar wall trace and initial condition use the original nominal mesh spacing. The physical wall mask uses the unregularized vertex distance. The floor affects both the inherited closure and the matching reaction split through their shared geometry array. The report includes the clipped fluid quadrature count and area fraction. The run retains existing convergence, positivity, scalar balance and tangent checks.

From the staged package parent submit `sbatch project4_fem/project43_wall_floor_sol.sbatch`. Results appear in `p43-wall-floor-<job_id>`.

## Decision after the matched Sol run (2026-11-24)

Array 2349247 used the same 0.5 mm mesh, Re = 50, solver revision and pseudo-transient continuation for both SST runs. The steady scaled residual tolerance was 1e-11.

| Task | Floor fraction | Final scaled residual | Interpretation |
|---|---:|---:|---|
| 2349247_1 | 0 | 2.61e-13 | Converged working baseline. Side-port share 4.57918%, static pressure loss 0.164719 Pa. |
| 2349247_2 | 0.5 | 1.64e-11 | Numerical near miss. Reached the minimum pseudo-time step, reported convergence false, and exited with code 2. |

Retain the original tolerance and failed status. The residual measures discrete equation imbalance. Its magnitude alone establishes neither the error in wall shear nor the accuracy of the split or pumping power. The stopping message does not identify which equation or constraint limited the update.

If an engineering comparison later needs the floor case, specify output stability tolerances before inspecting the last accepted iterates. Check wall shear, outlet split, pumping power, mass and scalar balances, interior k and omega, and the largest residual block. Any engineering acceptance applies only to the quantities checked; the recorded solver result remains a failure. Turbulent wall-shear verification against the flat-plate reference remains outstanding.

Defer further floor tuning. The next design stage minimizes laminar pumping power at Re = 50 from feasible Sol design 2349442, with fixed ports and inlet flow, side share 49.9–50.1%, and fluid area at most 83 mm².
