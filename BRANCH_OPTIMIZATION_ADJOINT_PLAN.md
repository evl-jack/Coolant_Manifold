# Coolant Manifold Branch Optimization and Adjoint Plan

**Draft for review**

Use a straight main bore with one radius shared by the inlet, pipe, and main outlet. Give the optimizer a separate side outlet radius and control over the side branch and its junction. Run the 2D search at a fixed inlet flow rate, then carry selected profiles to 3D OpenFOAM. Set the side/main target to **0.25** for the feasibility study and compare pressure loss and branch area relative to a fixed baseline.

Use **ASU Sol for all future 2D and 3D simulation and optimization runs**, including parameter sweeps and mesh convergence studies.

## Starting point and scope

The established workflow has already taken a 2D design through 3D geometry, meshing, and OpenFOAM RANS. The selected design reached **0.4005 side/main** on the 0.25 mm explicit-wall 2D grid. Its coarse 3D transient run averaged **0.2109 side/main** and **7.83 kPa** over 10–20 ms. The mean pressure drop increased from **7.67 to 7.99 kPa** across the two halves of that window. These are results for the earlier design with a downstream main-pipe pinch. The earlier 2D sensitivities showed that this pinch was the strongest local control on side flow, so reaching 0.25 with a straight main bore is the first question for the new design. [2D explicit-wall checks](rans_solver/results/shape_wall_checks.csv) · [Earlier 2D sensitivities](rans_solver/results/shape_gradients.csv) · [Sol transient summary](output/sol-urans-2061743/summary.json)

## Geometry controls

| Control | Proposed treatment |
|---|---|
| Shared main radius, Rₘ | One value per candidate, constant over the full 20 mm main passage. The main inlet and outlet use this same radius, with the bore open at both ends. Use 1.75 mm as the baseline radius because it matches the earlier main openings. |
| Side outlet radius, Rₛ | Independent design variable. Use 0.625 mm as the baseline radius, matching the earlier side opening. Keep the outlet plane and branch axis fixed during the first sweep. |
| Branch profile and junction | Start with a few smooth profile controls between the junction and outlet. Tie the branch root location to the main wall as Rₘ changes. The branch intersects the main at the junction, and the main bore remains open everywhere along its length. |
| Passage and manufacturing limits | Use a connected branch with one opening at each cross-section, symmetry about the branch axis, and bounds on minimum width and junction shape. Set the physical bounds during review. |

The shape search starts with a smooth half-width profile in the 2D plane. In 3D, revolve the main and branch profiles about their respective centerlines and join the fluid volumes. Use the same dimensions and profile definition in the 2D solver and the 3D builder. Additional branch controls can be added after the small profile model gives reliable sensitivities. [2D mesh](mesh.py) · [Shape search](rans_solver/shape_search.py) · [3D geometry builder](sol_rans_3d/build_mesh.py)

## Flow target and relative metrics

| Quantity | How it guides the design |
|---|---|
| Side/main flow ratio | Aim for 0.25 with net outward flow at both outlets. At mass balance this is 20% of inlet flow going to the side branch and 80% to the main outlet. |
| Relative pressure loss | Divide inlet-to-outlet static pressure drop by the pressure drop of the fixed baseline in the same model. Use the same pressure sampling method and inlet flow rate for every comparison. |
| Relative branch area | Divide branch fluid area outside the main bore by the baseline branch area. Track total 2D fluid area separately so we can see the cost of increasing the shared main radius. |
| Relative 3D volume | Measure fluid volume from the selected 3D geometry and compare it with the 3D baseline volume. Use the 3D flow and pressure results for final comparisons. |

Use the current 2D outlet pressure conditions of zero gauge pressure at both outlets for the feasibility study. Hold the 2D inlet flow per unit depth fixed. As the shared radius changes, recompute uniform inlet speed from the open inlet width. The proposed starting flow is **0.007 m²/s per unit depth**, corresponding to 2 m/s through the 3.5 mm baseline opening. For 3D, hold a separate volumetric inlet flow fixed across cases. A 3.5 mm round baseline inlet at 2 m/s gives approximately **1.155 L/min**.

Keep the existing signed mismatch, `e = (Q_side − 0.25 Q_main) / Q_in`, as the first search objective. It stays well defined near low main-outlet flow. Report the side/main ratio and outlet signs alongside it. The current acceptance band is `|e| ≤ 0.001`. Treat that as a starting proposal and agree on the tolerance for the actual 3D design at review.

Screen for the split first. Among candidates that meet the split, show the pressure-loss and branch-area tradeoff and retain representative designs for 3D. Begin with pressure and area as tradeoffs. A documented hardware limit can become a hard constraint. If the tested family falls short of the target, report the best feasible split and the limiting radius or branch constraint.

## Discrete adjoint and gradient checks

Write the converged 2D equations as a discrete residual R(w, d), where w contains the two velocity components, pressure, and the SST variables k and ω. The vector d contains the shared main radius, side outlet radius and branch profile controls. Solve a transpose system for each flow response J.

$$
R(w,d)=0,\qquad
R_w^\mathsf{T}\lambda=J_w^\mathsf{T},\qquad
\frac{\mathrm{d}J}{\mathrm{d}\mathbf d}=J_d-\lambda^\mathsf{T}R_d
$$

Use one adjoint for the signed flow mismatch and one for relative pressure loss. Where the main outlet flow is positive, recover the side/main ratio as `(0.25 + e) / (1 − e)` from the signed mismatch and fixed-flow conservation. Differentiate that expression to export its sensitivity. Use the analytic geometry derivative for branch area, including the dependence of the branch length on Rₘ. Save both the dimensional responses and the normalized gradients. [Giles and Pierce, *An Introduction to the Adjoint Approach to Design*](https://people.maths.ox.ac.uk/gilesm/files/ftc00.pdf)

Keep the background grid and residual layout fixed during optimization. Use differentiable cell coverage and fractional port openings. The inlet-radius derivative needs the changes in both open face area and inlet speed at fixed flow. Include the effects of wall distance, Brinkman resistance, face fluxes, and the SST closure in the residual derivatives used for the full RANS search. Use explicitly meshed walls to verify the candidate profiles after the search.

Check the laminar adjoint first, then repeat with the coupled RANS equations and the same geometry parameterization. Compare adjoint predictions against centered finite differences at several step sizes, with an inward difference at active bounds. Include a direction that changes more than one control. Tighten the flow solve until numerical noise is small relative to the perturbation response. Inspect both sign and magnitude. Near a vanishing sensitivity, compare absolute error in scaled variables. Report controls that cross a mesh or limiter transition so the optimizer can take a smaller step or use a checked local derivative. [NASA discrete-adjoint verification study](https://ntrs.nasa.gov/citations/19990113067)

## Work sequence and review evidence

1. **Build the shared geometry and baseline.** Define a single geometry record. Make the 2D and 3D builders consume it. Check the bore diameter at multiple stations, opening dimensions, connectivity and reconstruction of the fixed reference design.
2. **Screen the feasible geometry range in 2D.** Vary the shared radius and branch outlet radius and sample a few branch profiles. Record split, pressure-loss ratio and branch-area ratio. This establishes how close the straight-bore family can get to 0.25 before expanding the design space.
3. **Implement and verify the adjoint.** Check laminar flow, then coupled SST. Save gradient checks for each radius and branch control before enabling adjoint-driven search.
4. **Optimize and verify the shortlist.** Use the adjoint for the flow responses and analytical area derivatives. Re-solve selected candidates with explicit no-slip walls on two 2D grid resolutions. Check flow-split sensitivity to resolution before assigning a feasible label.
5. **Build and simulate selected 3D candidates on Sol.** Check mesh quality at each resolution. Compare the outlet flows, pressure drop, backflow and volume against the matching 3D baseline. If the run remains unsteady, compare successive averaging windows before reporting a final mean. Use the prior Sol case as a workflow reference. [Sol workflow](sol_rans_3d/README.md)

Save `candidate_results.csv` with geometry, baseline identifiers, all responses, convergence and feasibility for each evaluated design. Save `sensitivities.csv` with one row per design variable per candidate, including its value and units, the signed derivatives of flow ratio, relative pressure loss and relative branch area, adjoint status, and available finite-difference checks. Keep the geometry record with the shortlist so the dimensions used in 3D can be traced back to the exact 2D candidate.

## Decisions for review

- Confirm the allowable main and side radii, minimum branch width and junction shape constraints. The current opening dimensions provide baseline values.
- Agree on the flow-split acceptance band for 2D screening and for 3D verification. Confirm that both outlets should remain at equal prescribed pressure for the reference test.
- Confirm the inlet flow references and how to choose a candidate from the pressure-loss and branch-area tradeoff. Keep both values visible until that choice is made.

