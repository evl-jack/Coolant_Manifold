# Project 4 feasibility recovery and decision handoff

**Status on November 20, 2026.** This report evaluates the current fillet design and documents a proposed transition study. The **20% side-share target**, **19.5%–20.5% acceptance band**, existing **1.00 mm minimum passage width**, fixed port locations and openings, operating conditions, and current constraint definitions govern the approved case. The proposed transition and revised area definition require approval before an expanded FEM mesh or optimization is run. [Approved formulation and prior results](FORMULATION.md) · [Master project handoff](../COOLANT_MANIFOLD_MASTER_HANDOFF.md)

## Proposed transition experiment

Use two mirrored strips along the branch walls from **y = 6.500 to 8.000 mm**, centered about **x = 10.000 mm**. The original 20 fillet controls keep their existing support from **y = 5.000 to 6.500 mm**. At y = 6.500 mm the fixed throat spans **x = 9.500 to 10.500 mm**. The untouched terminal guide begins at y = 8.000 mm with edges at **x = 9.375 and 10.625 mm**. It extends to the side outlet at y = 11 mm. The proposed strips affect neither the straight main bore nor its exit guides.

For the new interval, set `t = (y − 6.5)/1.5` and define the existing branch half-width `h₀(t) = 0.5 + 0.125(3t² − 2t³)` mm. The proposed half-width is

\[
h(t)=h_0(t)+A(t)\sum_{i=1}^{6} B_i(t)(2d_i-1),
\qquad A(t)=16t^2(1-t)^2\,[h_0(t)-0.5].
\]

The six new independent controls lie in `[0,1]` and use cubic B-splines on the open knot vector `[0,0,0,0,1/3,2/3,1,1,1,1]`. Their Greville positions are **6.5000, 6.6667, 7.0000, 7.5000, 7.8333, and 8.0000 mm**. Apply the same coefficients to both walls by reflection. The endpoints in the list identify the first and last basis functions. Multiplication by `A(t)` fixes the actual wall endpoints and tangents. The combined design would have **26 independent controls**. The verification ladder must confirm the rank of all 26 on the new mesh.

The physical envelope on either wall is exactly `h₀(t) − A(t) ≤ |x − 10| ≤ h₀(t) + A(t)`. The largest half-width excursion is **0.07511 mm**, with overall physical x-extents approximately **9.3328–10.6672 mm**. For a later diffuse FEM embedding, the proposed support extends to `h₀(t) ± 1.14 A(t)` to include a cubic density transition of half-width `0.14 A(t)`. Its overall x-extents are about **9.3231–10.6769 mm**. Since the B-splines form a nonnegative partition of unity and `A ≤ h₀ − 0.5`, the minimum physical width is at least **1.00 mm**. Each cross-section remains a single interval and attaches to the original connected passage. The full T would retain its directional inlet and solve without imposing flow symmetry. [Proposed and existing envelopes](output/sol-2083406/design-envelopes.png)

If this experiment is approved, first test common narrowing and widening at several amplitudes, upper-versus-lower curvature changes, and a compensating fillet-plus-transition mode. Run Stokes and Re 50 at the same inlet flow, traction conditions, and held dimensions. Measure directional derivatives at the starting design. Use actual simulated response ranges and geometry checks when assessing feasibility.

## Area accounting for the proposed region

Keep the physical intent of the existing cap, which is to preserve **net editable fluid area**. Recompute the baseline for the expanded region on the matched mesh. Over the union of the fillet and transition density strips, require `∫(γₕ,baseline − γₕ,design) dA ≤ 0`. In the reconstructed geometry, report the corresponding signed change `2∫[h₀(y) − h(y)] dy` from **y = 5 to 8 mm**. The two calculations both include the onset of the fillet at y = 5 mm. The baseline-relative formulation accounts for the added support area and permits local narrowing accompanied by compensating widening. The old **0.396621 mm²** solid cap continues to describe the original strips.

## Stage A results on Sol

Job **2083406** reproduced all **16 archived FEM states**. For each model, this included the seven archived samples and the constrained range probe. The largest absolute difference across the four archived scalar responses, each in its reported units, was **7.1×10⁻¹²**. All replayed states met the solver and mass-balance checks. The table includes the original-resolution controls and an additional two-level resolution check at baseline. Pressure and dissipation are nondimensional. Leakage is RMS speed at quadrature points with `γₕ < 0.05` as a percentage of mean inlet speed. [All replayed cases](output/sol-2083406/replayed-cases.csv) · [Audit results](output/sol-2083406/audit.json)

| Case | Side share | Pressure loss | Dissipation | Solid area (mm²) | Minimum width (mm) | Solid RMS speed | Converged |
|---|---:|---:|---:|---:|---:|---:|---|
| Stokes baseline | 27.648% | 107.9078 | 107.5663 | 0.396621 | 1.000 | 0.654% | Yes |
| Stokes range probe | 27.971% | 106.7311 | 106.3881 | 0.048999 | 1.000 | 0.367% | Yes |
| Stokes fine baseline | 27.629% | 108.0559 | 107.4801 | 0.396480 | 1.000 | 0.645% | Yes |
| Re 50 baseline | 25.101% | 90.7808 | 111.1367 | 0.396621 | 1.000 | 0.984% | Yes |
| Re 50 range probe | 24.313% | 92.2899 | 112.2720 | 0.395229 | 1.000 | 1.212% | Yes |
| Re 50 fine baseline | 25.089% | 90.7474 | 111.0933 | 0.396480 | 1.000 | 0.944% | Yes |

The Stokes probe's flow state converged while its pressure exceeded the limit. The Re 50 probe's pressure and area satisfy their limits and its side share exceeds the target band. The earlier boundary-fitted re-solve of that probe gave **24.0707%** on its fine mesh. [Boundary-fitted diagnostic](output/sol-2082990/explicit-probe.json) · [Frozen Re 50 FEM velocity](output/sol-2083406/baseline-velocity.png)

### What the geometry can influence

On the existing fillets, the largest baseline Re 50 side-share gradient was **control 8**, about **+0.2025 percentage points per unit coefficient**. Controls **9** and **7** followed, at about **+0.1852** and **+0.1723**. Their Greville positions lie in the lower fillet around **y = 5.16–5.30 mm**. In Stokes, controls **15, 14, and 16** were strongest, at about **+0.2910**, **+0.2799**, and **+0.2451 percentage points per unit coefficient**, around **y = 5.83–6.09 mm**. These are derivatives at the baseline. The five coordinated fillet perturbations actually gave **24.1205%–25.5272%** at Re 50 and **26.5751%–27.9709%** in Stokes. The all-narrowed points exceed the original area allowance.

For the proposed interval, a **friction-only screening calculation** estimates branch resistance proportional to `∫dy/w³` from y = 5 to 11 mm for a planar slit of width `w`. The estimate holds the effective main-path resistance and junction driving pressure relationship fixed and is calibrated to the **24.0707% explicit Re 50 probe**. It omits inertial turning losses, geometry-induced flow redistribution at the junction, and area compensation. At maximum coordinated narrowing, the proposed interval adds approximately **6.37%** branch resistance and moves that estimate to **22.96%** side share. The model needs about **22.94%** more branch resistance to reach **20.5%**, or **26.81%** to reach **20.0%**. Even the idealized friction comparison of holding the entire 1.5 mm transition at 1.00 mm gives **22.01%**. That comparison disregards the smooth joins and area constraint. These figures are analytical scenarios for judging the pilot's value. They provide neither an expanded-geometry FEM result nor a bound on the full PDE. Full coordinated narrowing would remove about **0.1000 mm²** of fluid area and would require compensation under the proposed area rule. [Measured response and analytical amplitude plot](output/sol-2083406/side-flow-authority.png)

## Second hydraulic constraint and shared pressure cap

On the fixed-flow FEM model, the boundary reaction work balances dissipation plus the assembled advective work, with the weak discrete pressure work near roundoff. The largest work-balance residual across the replay was under **3.1×10⁻¹²**. In the Stokes baseline, reaction work and dissipation both equal **107.5663** to the shown precision. In Re 50, reaction work is **90.3349**, dissipation **111.1367**, and advective work **−20.8019**. The current pressure loss uses an area-mean inlet static pressure and a signed-flow-weighted static outlet pressure. Its numerical values are **107.9078** and **90.7808**. The energy balance explains why pressure and dissipation move closely in this fixed-flow setup. The specific correlation depends on the geometry, pressure averaging, and changes in advective work.

| Criterion | Current static pressure loss | Main-outlet profile nonuniformity `Jᵤ` |
|---|---|---|
| Physical use | A pressure allowance at specified reference planes can represent a real head constraint when tied to hardware and the operating point. | A profile limit can protect a sensitive downstream component when there is a requirement at the measurement plane. |
| Relative design-gradient norm at baseline | **0.00783** Stokes and **0.01012** Re 50 per unit control-vector norm. | About **3.0×10⁻¹²** Stokes and **2.54×10⁻⁵** Re 50. |
| Correlation with dissipation across seven samples | **0.999997** Stokes and **0.9947** Re 50. | Stokes variation was about **10⁻¹²**, making the sample correlation uninformative. Re 50 correlation was **0.899** across its small range. |
| Conditioning and resolution | Baseline pressure shifted by **0.148** in Stokes and **0.033** in Re 50 under refinement. Archived pressure gradients passed the manifold checks. | At Re 50, the sample range was **1.61×10⁻⁵** while refinement changed the baseline by **8.15×10⁻⁴**, about **51 times** more. In Stokes, refinement changed it by **0.00577**, much larger than the observed design variation. |
| Cost | Boundary integration and one adjoint right-hand side. | Exact P2 outlet integration and one adjoint right-hand side. Re 50 directional checks at two step sizes agreed within **2.3×10⁻¹⁰** absolute. The current outlet location offers little verified sensitivity. |
| Common feasibility | All sampled Stokes cases exceeded the current **95.3198** cap. The best observed loss **106.7311** was about **12%** higher. | No engineering threshold has been specified, so feasibility remains open in both models. |

**Recommendation.** Retain pressure as the second hydraulic response and prefer a **physically specified pressure allowance at fixed reference planes** for the constraint. The current main-outlet nonuniformity is a poor replacement because mesh changes overwhelm its variation with the design. Set a common pressure limit only when hardware requirements and the pressure definition support one at the fixed inflow. The present **95.3198** cap comes from the Re 50 baseline and has no independent hardware basis. A common cap is defensible if it represents the same device requirement for both models, even if one model then predicts infeasibility. Stokes here is the same reduced setup with convection removed. Its larger pressure is a result to report, rather than a reason to select a new cap. Approval and a physical requirement are needed before using a different cap or response. [Force balance, mesh comparison, and metric derivatives](output/sol-2083406/audit.json)

## Cost and approval record

Job **2083406** took **69 seconds** on one Sol CPU, including the two-grid baseline checks, candidate metric adjoints, and figure generation. Eight archived cases took **1.78 s** in Stokes and **8.76 s** at Re 50 on **2,116 triangles**. The refined baseline on **8,464 triangles** took **1.73 s** and **7.83 s**, respectively. A proposed initial mesh could raise transition resolution from six to eighteen rows, adding about **384 triangles** for roughly **2,500** total. That supports a planning allowance of **1–3 CPU minutes** for a few dozen baseline-resolution cases and gradient checks in both models. Benchmark the expanded implementation before estimating full optimization runtime. Each adjoint response needs one transpose right-hand side, and the two flow-band inequalities share one response. A centered finite-difference gradient over 26 controls would require **52 primal solves per design**.

| Item | Original requirement and issue | Proposed action and expected impact | Approval status |
|---|---|---|---|
| Editable geometry | Twenty controls in the fillets from y = 5 to 6.5 mm. The best explicit probe remains 3.57 points above the allowed side share. | Add six controls from y = 6.5 to 8 mm, preserving the fixed geometry and 1.00 mm width. Test the estimated modest gain before deciding on a broader region. | **Proposed.** Analytical and geometric assessment complete. Expanded FEM requires approval. |
| Area | Baseline-equivalent solid cap of 0.396621 mm² in the original strips. | Use the matched expanded baseline and zero net fluid-area loss over the union. This keeps the meaning of the budget while enabling compensating shape changes. | **Proposed.** Approval required before applying it. |
| Pressure | Static-pressure cap of 95.3198 based on Re 50. All sampled Stokes points exceed it. | Retain it for provenance and screening. Request the operating pressure requirement and its measurement definition before adopting a common optimization cap. | **Assessment complete.** A cap change requires approval. |
| Second hydraulic requirement | Pressure loss is the existing constraint. The outlet uniformity candidate has poor design signal at the current outlet plane. | Prefer a justified pressure allowance. Consider a different physical response only when its measurement location, engineering need, and resolved design sensitivity are available. | **Recommendation only.** Substitution requires approval. |
| Target, tolerance, width, and physics | 20% side share with ±0.5 percentage points, 1.00 mm width, Stokes and Re 50 at the same conditions. | Carry these through the proposed pilot and preserve all fixed ports and collars. Reuse the current gradient and convergence acceptance criteria. | **Existing requirements.** Any revision requires approval. |

The next executable stage after approval is a small matched Stokes and Re 50 transition study, followed by map and adjoint checks for all 26 controls. Use measured responses to decide whether the full optimization, KKT analysis, and optimized-design presentation are warranted.

## Reproduce the evidence

Use a Sol compute allocation with the scientific Python environment. Replay inputs live in [the archived fillet study](output/sol-2082990/search-re-50.json) and [the frozen FEM reference](output/sol-2082569/junction-reference.json). The audit adds diagnostic responses and an analytical transition estimate to those conditions. [Audit source](recovery_audit.py) · [Plot source](plot_recovery_audit.py)

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -m project4_fem.recovery_audit --archive project4_fem/output/sol-2082990 \
  --output project4_fem/output/recovery-replay
python -m project4_fem.plot_recovery_audit --output project4_fem/output/recovery-replay \
  --baseline project4_fem/output/sol-2082569/reference-grid-2.npz
```
