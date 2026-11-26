# Fillet optimization decision record

## Frozen case and design map

Use the verified P2–P1 solver, with the existing parabolic inlet, prescribed unit nondimensional inlet flow per depth, no-slip fixed walls, and natural zero traction at both outlets. Use Re 50 for the nonlinear design and switch off convection for the Stokes comparison. The main bore, branch throat, straight runs, and all three openings use the frozen Sol geometry. [Reference geometry and flow](output/sol-2082569/junction-reference.json)

The editable region consists of the two strips along the circular fillet between the main bore wall at **y = 5 mm** and the throat at **y = 6.5 mm**. Let `s` be the normalized angle around the original quarter-circle, increasing from the main wall to the throat. The design has **20 independent cubic B-spline coefficients** in `[0,1]`. The same coefficients act on both fillets by reflection about **x = 10 mm**. The flow is solved on the entire T.

The branch half-width is `w = w₀ + a(y) Σ Bᵢ(s)(2γᵢ − 1)`. The maximum raw half-width change is `0.45 sin²(πs)` mm. Combine this with the clearance `w₀ − w_throat` by their harmonic product `a = raw × clearance / (raw + clearance)`, with zero at the joins. The displacement and its angular derivative vanish at both joins, and the clear width stays at least as large as the fixed throat. A compact cubic transition of half-width `δ = 0.14a` gives `γₕ = 3z² − 2z³`, where `z = clip(0.5 + (w − |x − 10 mm|)/(2δ), 0, 1)`. The fixed computational mesh covers the entire outward strip. Beyond the strip the supported domain is fixed and has zero drag. The analytical contour for reconstruction is `|x − 10 mm| = w`.

## Responses and limits

| Item | Definition and scale |
|---|---|
| Objective | Minimize the sum of viscous and Brinkman dissipation `∫ [2με:ε + α|u|²] dΩ`, divided by its Re 50 baseline value. |
| Flow distribution | Keep `fₛ = Qₛ/Qᵢₙ` within `0.20 ± 0.005`. The default side/main target `0.25` is configurable. Scale the two inequalities by `0.005`. |
| Pressure loss | Use the area-mean inlet pressure minus the signed-flow-weighted mean outlet pressure across both fixed outlet planes. The pressure trace and normal velocity come from the same discrete state. Scale by the cap in the constraint `Δp / Δpmax − 1 ≤ 0`. |
| Pressure cap | Set `Δpmax` to **105% of the Re 50 Brinkman baseline** on the design mesh. This course screening cap permits a small extra pressure cost to meet the flow target. A hardware head budget would replace this provisional value. Use the same cap for Stokes and Re 50. |
| Area | Integrate `1 − γₕ` over both complete editable strips and require it to be at or below the baseline solid area. Scale by that positive baseline area. This permits redistribution while preserving net fluid area in the edit region. |

Use `α = 20,000 (1 − γₕ)²` in nondimensional units. Record a baseline penalty check at `5,000`, `20,000`, and `80,000` and measure speed inside strongly penalized material. This check addresses leakage and resolution sensitivity. Possible changes to the penalty based on these data require an explicit record. Keep the fixed state mesh through optimization and re-solve shortlisted contours on boundary-fitted grids.

For each response, verify the adjoint with centered differences in two random directions and a collective fillet displacement. Use steps `0.04`, `0.02`, `0.01`, and `0.005` from the interior baseline. Record absolute and relative errors and the Taylor remainder. Require directional agreement within `max(10⁻⁷, 0.002 |D|)` at the finest step and observed second-order remainder over a usable interval. Require the free-degree Newton residual to meet the solver tolerance and relative port mass imbalance below `10⁻⁹`. Report feasibility using the declared band and limits with scaled numerical tolerance `10⁻⁶`. Use the actual scaled constraints in the KKT calculation.

## Change record

| Original proposal | Issue and selected implementation | Consequence and status |
|---|---|---|
| Local smooth density mapping with at least 20 controls and justified mirroring | The support and the mirror plane need numerical definitions. Use the angle-based half-width map and paired strips above. | The flow passages remain connected by construction, with fixed joins and a fixed throat. Numerical map and hydraulic checks are required. **Selected for implementation.** |
| Pressure loss first, with an engineering rationale for its cap | Equal outlet traction still allows a nonzero outlet pressure trace. Use the fixed-section average above with a baseline-relative cap. | Report correlation and gradient alignment with dissipation. Pressure remains the second PDE-dependent requirement, even if the cap is inactive. **Selected for screening.** |
| A solid-area limit or minimum fluid-area requirement | On this fixed two-strip region, solid and fluid-equivalent areas sum to the region area. | A cap at baseline solid area implements no net loss of editable fluid area across the full fillet. **Selected for screening.** |
| Candidate Brinkman polynomial and optional scaling | The penalty, exponents, and scales were unspecified. Use the constants and baseline normalization above. | The penalty study measures leakage and sensitivity before interpreting the physical contour. **Selected for screening.** |
| Verification steps from 0.04 through 0.005 | The first Re 50 random direction passed the centered-difference tolerance but its side-share Taylor slopes were 3.05, 3.42, and 4.48 over that interval. Add `0.0025`, `0.00125`, and `0.000625` to check the asymptotic behavior after the apparent cancellation. | The same directional and Taylor acceptance limits apply. **Implemented after Sol job 2082955.** |

The optimization target, tolerance, Reynolds number, pressure-constraint type, and boundary conditions retain the proposal's intent. A finite feasibility search can establish an observed feasible point. If every result misses the target, preserve the problem definition and report the search evidence before considering a revised engineering question.

## Sol findings

Job **2082990** verified the design map and all four response gradients in Stokes and Re 50. The design tangent had rank 20, the reflected density fields agreed, and the minimum sampled passage width was **1.00 mm**. The largest relative directional error at step `0.000625` was `8.64e-6`. The extended Taylor sequence included an interval of order 1.84 for the previously unresolved side-share direction. [Map checks](output/sol-2082990/map.json) · [Re 50 gradients](output/sol-2082990/gradients-re-50.json) · [Stokes gradients](output/sol-2082990/gradients-re-0.json)

The realized caps were **95.3198** in nondimensional pressure units and **0.396621 mm²** of editable solid. The Re 50 samples spanned **24.1205%–25.5272%** side share. The **24.1205%** case exceeded the area cap. The constrained range probe reached **24.3133%** with pressure loss **92.2899** and solid area **0.395229 mm²**, then reached its iteration limit. The best Stokes pressure observed in the initial samples was **106.7311**, **12.0%** above the common cap. The search found zero cases in the **19.5%–20.5%** target band that met all limits. These are local feasibility results from a finite search. [Limits](output/sol-2082990/limits.json) · [Re 50 search](output/sol-2082990/search-re-50.json) · [Stokes search](output/sol-2082990/search-re-0.json)

Pressure and dissipation had sample correlations of **0.9947** at Re 50 and **0.999997** in Stokes. The baseline gradient cosine was **0.9952** at Re 50 and nearly one in Stokes. Pressure served as a separate head-budget requirement and excluded the sampled Stokes cases. These data provide little evidence of an independent hydraulic tradeoff across this design region. Any replacement should have a physical specification and user approval before another optimization. [Re 50 search](output/sol-2082990/search-re-50.json) · [Stokes search](output/sol-2082990/search-re-0.json)

The baseline Re 50 side share was **25.2133%**, **25.1011%**, and **25.0102%** at the three trial penalties. RMS speed in nominal solid fell from **2.68%** to **0.34%** of the inlet mean speed over the same range. The selected penalty stayed at **20,000** for the feasibility experiment. [Penalty study](output/sol-2082990/penalty.json)

Job **2083052** re-solved the baseline and the Re 50 range-probe contour on boundary-fitted meshes of 2,116 and 8,464 triangles. The fine baseline side share was **24.8960%**. The fine probe gave **24.0707%** side share, pressure loss **92.7530**, and geometric solid area **0.394982 mm²**. The probe's side share changed by **0.0215 percentage points** between the fitted meshes. It remained **3.5707 percentage points** above the upper flow limit. The maximum realized half-width excursion allowed by the map is **0.2455 mm**. [Explicit-wall probe](output/sol-2082990/explicit-probe.json)

**Decision.** The verified design map can shift flow, and the tested fillet-only geometry fell short of the 20% target. The common pressure cap also excluded the sampled Stokes designs. Optimization and KKT reporting require a feasible engineering specification. If the 20% target is mandatory, obtain approval to study a larger editable region while retaining the fixed ports, and review the shared pressure requirement. If the fillet-only region is mandatory, establish a reachable common flow band and a justified second hydraulic requirement before optimizing either model.

To reproduce the gated checks and the diagnostic in a Sol compute allocation, run from the project root:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -m project4_fem.fillet_study --output project4_fem/output/fillet-study
python -m project4_fem.fillet_explicit --study project4_fem/output/fillet-study \
  --output project4_fem/output/fillet-study/explicit-probe.json
```
