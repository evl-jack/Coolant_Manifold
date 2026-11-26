# Coolant manifold common case specification

Frozen for the bounded next stage, before independent execution. The reference records are Sol jobs **2083835** and **2084088**. The exact baseline and screened-candidate coefficient arrays are in the [matched result manifest](output/sol-2084088/matched-verification.json). Preserve the source archives. Use independent output directories and record the copies of the input manifests used by each job.

## Geometry and state

- Coordinates and all shape definitions are in **mm**. The planar inlet is x = −1 at y = 3–5, the main outlet is x = 21 at y = 3–5, and the side outlet is y = 11 at x = 9.375–10.625. The main bore is 2.00 mm high and the side branch is centered at x = 10. The saved 26-control fitted walls follow the analytical midpoint contour. Curves may be approximated by straight mesh edges, with equivalent refined levels. Keep the fixed mouths, bores, and openings for every paired comparison.
- Use water ρ = 998.2 kg/m³ and μ = 1.002×10⁻³ Pa·s to dimensionalize the laminar OpenFOAM reference. Define H = 2.00 mm and U = Re μ/(ρH). At Re = 50, U is approximately 0.025095 m/s. Use `6 U η(1−η)` along +x at the inlet, where η = (y−3 mm)/H. Its mean is U and prescribed flow per depth is UH. Use no-slip solid walls. An OpenFOAM one-cell extrusion of 1 mm with front/back `empty` is a numerical representation of a planar domain. That thickness supplies only the flow and integration scale and is separate from hardware depth.
- FEM uses natural zero outlet traction for the symmetric Newtonian stress, with no prescribed outlet velocity. Record the equations, traction approximation, and any imposed outlet pressure/velocity gradient of the OpenFOAM counterpart. In both codes use signed outward main/side flux and a positive inlet flux. Check clipping or reverse flow at the outlets.
- Compare velocity normalized by U at common interior sections x = 0 and 20 mm and branch y = 10.5 mm for the archived pair. Sample both components where feasible. The integral response planes remain the three port boundaries. Normalize OpenFOAM kinematic pressure using νU/H and dynamic pressure using μU/H. Compare area-mean inlet static pressure minus combined signed-flow-weighted outlet static pressure. Normalize dissipation per unit depth by μU² and compare `∫2μ S:S dA` only when evaluated consistently.
- Record side share `Qs/Qin`, side/main `Qs/Qm`, relative flux imbalance `|Qin−Qm−Qs|/Qin` and physical and integrated area. A positive candidate benefit is **baseline side fraction minus candidate side fraction** in percentage points within one method, mesh level and operating condition. Analytic area change is the signed coolant area difference over the complete changed support relative to the **26-control baseline**. Positive is added fluid. A positive inequality residual is a violation. Bounds stay **19.5%–20.5%**, **1.00 mm minimum width**, **95.31983921402909** provisional nondimensional pressure cap, and no net area loss for this pilot. Geometry, pressure, area and flow flags stay separate.
- Mesh naming is `L1`, `L2`, `L3` for three systematic levels. Record actual counts and representative wall spacings. Re 50 is laminar. Stokes data are eligible for direct comparison only if convection is actually disabled. Use fitted-wall, Brinkman and OpenFOAM result labels explicitly. Do not infer turbulent operating conditions from the laminar case, the numerical extrusion depth, or the old 2 m/s setup.

## Predeclared reporting gates

FEM retains its established verification checks. For OpenFOAM archive the native residual history, port balance, mesh checks, and stability of flow split and pressure. Aim for relative mass imbalance ≤10⁻⁶ and stable last-window side share within 0.005 percentage points, pressure within 0.1%, and mesh level 2-to-3 differences within 0.10 points and 0.5% respectively. For cross-code corroboration, use thresholds of 0.20 points for the two fine side fractions, 0.10 points for the paired gain, 2% for static pressure, and 3% of U for velocity-profile RMS differences where sampled. Flag misses and quantify the outlet model discrepancy. These numerical thresholds do not alter the hydraulic design constraints.

The turbulent reference uses fitted walls. State its independent hydraulic diameter, roughness and inlet development assumptions. Record actual hardware inputs or mark assumed benchmark inputs as such. Audit the installed OpenFOAM SST variant against the planned SST-2003m equations, coefficients and pressure convention. Apply the RANS residual, Jacobian and derivative gates in the [prior decision handoff](GEOMETRY_AND_RANS_DECISION_HANDOFF.md). Validation of a benchmark or one residual layer does not establish whole-manifold or full adjoint validation.

## Ownership

| Owner | Exclusive source area and Sol run prefix |
| --- | --- |
| Lead | Common specification and integrated decision report under `project4_fem`. Turbulent reference definition and case under `project4_fem/next_turbulent_reference`. |
| A — OpenFOAM | `project4_fem/next_openfoam` and scratch `p4-laminar-foam.*`. Owns laminar cases, mesh comparison and cross-checks with B and C. |
| B — geometry | `project4_fem/next_geometry` and scratch `p4-metering.*`. Owns preregistered DOE, as-built geometry checks, and fitted-wall screening. |
| D — RANS adjoint | `project4_fem/next_rans_adjoint` and scratch `p4-rans-adjoint.*`. Owns fitted RANS foundations and derivative verification. |
| E — HPC | `project4_fem/next_hpc` and scratch `p4-hpc.*`. Owns environment, manifests and Slurm cost accounting. Coordinate MPI tests with A using a ready case. |

Existing FEM, shape and solver interfaces are shared read-only for this stage. Use unique scratch roots and account for failed runs. Record input hashes, commit where accessible, environment, resource requests, elapsed time and peak memory. Maintain a complete record of gates that pass, fail, or lack evidence.
