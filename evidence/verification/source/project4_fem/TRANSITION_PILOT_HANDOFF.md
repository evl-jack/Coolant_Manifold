# Project 4 handoff — 26-control FEM feasibility pilot

**Decision.** The six transition controls lowered side flow, but the tested designs did not demonstrate enough authority to justify a full optimization toward the **19.5%–20.5%** band. The closest screened Re 50 case delivered **24.2258%** of inlet flow to the side outlet. It met the pressure, area and width limits, and missed the upper flow limit by **3.7258 percentage points**. The pilot establishes a local design effect on the tested mesh. Feasibility over the full 26-dimensional design space remains open.

## Geometry and operating point

The pilot retains the 20 fillet controls over `y = 5.0–6.5` mm and adds six independent cubic B-spline coefficients over `y = 6.5–8.0` mm. Both walls reflect about `x = 10` mm. With $t=(y-6.5)/1.5$, the original transition half-width is $h_0=0.5+0.125t^2(3-2t)$ mm. The approved adjustment is $A(t)\sum_i B_i(t)(2d_i-1)$, with $A=16t^2(1-t)^2(h_0-0.5)$, coefficients $d_i\in[0,1]$, and knots `[0,0,0,0,1/3,2/3,1,1,1,1]`. The peak wall excursion on either side is 0.07511 mm. The compact cubic density transition preserves the prior fillet law and uses the same $0.14A$ half-width in the new support.

The map passed before each PDE sweep. The density Jacobian had rank **26**, with smallest singular value **0.03216**. The two map derivative checks had maximum pointwise errors of **3.18×10⁻⁹** and **3.75×10⁻⁹**. Fixed-region and mirrored-density errors were zero on the checked meshes. Wall displacement at the joins was zero. The analytical envelope gives zero endpoint slopes, and the sampled one-sided tangent defect was below **8.9×10⁻⁶**. The cubic density ramp reaches zero slope at saturation, with densities in `[0,1]`. A continuous, single-interval branch section at every station joins the continuous main bore. The nonnegative spline basis and envelope preserve at least **1.00 mm** total width throughout the design bounds. Inlet and main outlet widths remain **2.00 mm** at `x = −1` and `x = 21` mm. The side outlet remains **1.25 mm** wide at `y = 11` mm, beyond the fixed guide starting at `y = 8` mm. [Map verification](output/sol-2083835/map-re-50.json) · [Geometry formulation and prior audit](RECOVERY_DECISION_HANDOFF.md)

Both models used the same 2,500-triangle Taylor–Hood FEM mesh, prescribed parabolic inlet with unit nondimensional flow per depth, natural zero traction at both outlets, and Brinkman coefficient 20,000. The pilot increased the number of transition rows from six to 18 to resolve the added strip. It uses the reduced 2D Stokes and Re 50 equations from Project 4. The target remains $Q_s/Q_{in}=0.20$, equivalent to $Q_s/Q_m=0.25$ under mass conservation.

## Measured hydraulic response

Each model used the same 13 designs, including four narrowing levels, two widening levels, area-neutral curvature redistribution, the archived Re 50 fillet probe, and four compensated combinations from the baseline and archived fillets. All **26 of 26** FEM states converged. The maximum nonlinear residual was **9.73×10⁻¹⁴** and the maximum relative mass imbalance was **1.38×10⁻¹³**. Every sampled design retained the 1.00 mm throat. The complete records include the 26 coefficients, flow ratio, relative mass imbalance, nonlinear residual, dissipation, leakage, both area measures, and signed constraint residuals for each case. [Complete results CSV](output/sol-2083835/pilot-cases.csv) · [Stokes case records](output/sol-2083835/cases-re-0.json) · [Re 50 case records](output/sol-2083835/cases-re-50.json)

| Model | Baseline side share | Measured range | Closest case meeting area and width | Pressure loss for closest case | In target band |
| --- | ---: | ---: | ---: | ---: | ---: |
| Stokes | 27.9587% | 26.8299%–28.4897% | 26.9241% | 108.6881, above cap | 0 of 13 |
| Re 50 | 25.3151% | 24.2258%–25.6693% | **24.2258%** | **92.0934**, below cap | 0 of 13 |

The closest area-compliant case in both models paired the archived Re 50 fillet shape with **full transition narrowing** and compensated by widening the fillet. At Re 50 it gave $Q_s/Q_m=0.319710$, dissipation **112.0329**, and nominal-solid RMS leakage of **1.079%** of mean inlet speed. Its side share fell **1.0893 percentage points** from the matched baseline and **0.2918 points** from the same-mesh archived fillet probe. From the baseline fillet alone, full narrowing lowered the Re 50 side share by **0.7921 points**. Restoring its area reduced that gain to **0.6234 points**. The area-neutral curvature mode changed Re 50 side share by only **−0.0335 points**.

The original 20-control Re 50 Brinkman probe gave **24.3133%** under its own area rule. The new best result is **0.0875 points** lower. The historical boundary-fitted re-solve gave **24.0707%** for the old probe. The expanded baseline differed from the old Brinkman baseline by **+0.2140 points** at Re 50 and **+0.3104** in Stokes. The expanded pilot used a single mesh and the Brinkman wall treatment, so the matched within-pilot change is the relevant evidence for added geometry authority. The physical gain over the original 20-control case remains unverified by a matched grid and wall-treatment comparison. [Historical audit and fitted-wall reference](RECOVERY_DECISION_HANDOFF.md) · [Historical Re 50 probe](output/sol-2082990/search-re-50.json) · [Historical fitted-wall case](output/sol-2082990/explicit-probe.json)

### Analytical comparison

The friction-only proxy scales branch resistance as $\int dy/w^3$ and is recalibrated to the matching FEM baseline in each model. At full coordinated narrowing, it predicts **26.7258%** side share for Stokes and **24.1597%** for Re 50. FEM measured **26.8299%** and **24.5230%** in those same cases. The proxy overpredicted the narrowing effect by **0.1040** and **0.3633** percentage points, respectively. The older **22.96%** estimate used the historical explicit-wall Re 50 probe as its calibration and omitted area compensation. It remains a friction estimate. The measured results capture junction effects and fluid redistribution that the proxy holds fixed. [Analytical comparison](output/sol-2083835/prediction-versus-fem.pdf)

## Area and pressure screens

Positive area change means added fluid area. The diffuse measure integrates `γ_design − γ_baseline` over the paired editable strips on the matched mesh. The physical measure is $2\int_5^8 [h(y)-h_0(y)]\,dy$, evaluated from the reconstructed wall. Each uses the same baseline and full `y = 5–8` mm support. The expanded baseline has **0.509994 mm²** of diffuse solid area within the strips. The historical **0.396621 mm²** cap remains tied to the original fillet strips. Nine designs per model preserved baseline area in both measures, and all four uncompensated narrowing levels lost area. Across the designs, the greatest discrepancy between the two area measures was **0.001019 mm²**.

| Design | Diffuse fluid-area change (mm²) | Reconstructed fluid-area change (mm²) | Diffuse minus reconstructed (mm²) |
| --- | ---: | ---: | ---: |
| Full narrowing | −0.100708 | −0.100000 | −0.000708 |
| Full widening | +0.099862 | +0.100000 | −0.000138 |
| Curvature redistribution | +0.000000 | +0.000100 | −0.000100 |
| Archived fillet with full narrowing and compensation | +0.000000 | +0.000820 | −0.000820 |

Pressure loss is the **area-mean static inlet pressure minus the signed-flow-weighted static pressure over both outlets**, at the same fixed planes and inlet flow. Dissipation combines viscous and Brinkman drag work. Both are nondimensional. Across the 13 cases, their Pearson correlation was **0.99999** in Stokes and **0.99395** at Re 50. The earlier energy audit identifies the contribution of advective work at Re 50 and explains why the static pressure metric and dissipation track one another in this fixed-flow setup. The provisional pressure cap of **95.3198** admitted all 13 Re 50 cases (**90.1031–92.0934**) and excluded all 13 Stokes cases (**107.3885–108.7945**). The smallest Stokes pressure loss exceeded the cap by **12.0686**. Keep pressure as the second hydraulic constraint candidate. The prior mesh audit found that changes in outlet-profile nonuniformity due to grid refinement swamped its design signal at the current plane. The physical pressure budget and its reference planes still require a hardware specification. [Energy and outlet metric assessment](RECOVERY_DECISION_HANDOFF.md)

## Adjoint and convergence verification

The discrete adjoint uses the 26-column density tangent and the existing PDE linearization. The side/main ratio follows from the mass-conserving side fraction through $f_s/(1-f_s)$. Both area responses have direct design derivatives. Nominal-solid leakage uses a density threshold $\gamma<0.05$ for reporting. Leakage, mass imbalance, and nonlinear residual are state diagnostics.

At the baseline, both models passed the same centered finite-difference and Taylor remainder tests for every response in two random directions, one all-control direction, and one transition-only direction. The seven step sizes were **0.04** through **0.000625**, halving each time. Each perturbed state met residual ≤ 10⁻¹¹ and relative mass imbalance < 10⁻⁹. The original acceptance checked the final centered error against `max(1e-7, 0.002·|adjoint directional derivative|)`. It also required a sampled Taylor order in `(1.7, 2.3)` or a final remainder below `1e-12`. The largest final error relative to the permitted tolerance was **0.57%** for Stokes and **1.50%** at Re 50. The largest Taylor remainder fell from **1.88×10⁻⁴ to 3.60×10⁻⁸** in Stokes and **2.75×10⁻⁴ to 5.32×10⁻⁸** at Re 50. The full step histories show small non-monotonic errors near the smallest steps.

| Response | Stokes maximum final centered error | Re 50 maximum final centered error |
| --- | ---: | ---: |
| Dissipation | 3.52×10⁻⁶ | 2.67×10⁻⁶ |
| Side fraction | 5.66×10⁻⁸ | 3.60×10⁻⁸ |
| Pressure loss | 3.49×10⁻⁶ | 4.30×10⁻⁶ |
| Diffuse solid area | 8.43×10⁻⁷ | 8.43×10⁻⁷ |
| Side/main ratio | 1.09×10⁻⁷ | 6.46×10⁻⁸ |
| Reconstructed fluid-area change | 1.34×10⁻¹⁴ | 1.34×10⁻¹⁴ |

[Stokes derivative records](output/sol-2083835/gradients-re-0.json) · [Re 50 derivative records](output/sol-2083835/gradients-re-50.json)

## Next formulation and decisions

The Re 50 shortfall is primarily limited by the sampled passage authority. Area compensation reduced the response of coordinated narrowing, and Re 50 pressure had headroom in every case. In Stokes the shared pressure cap adds a separate obstacle. Stop at this pilot. Before approving a run with an expanded design space, decide which additional geometric freedom is physically acceptable and whether Stokes is a screening model or must satisfy the same hardware cap as Re 50. The pilot's area rule also needs explicit approval for production use.

If a new region is approved, first verify its map, then compare a matched baseline and screened candidate on refined and boundary-fitted meshes. For a later constrained optimization, keep the flow band as a hard constraint, apply the approved baseline-relative area and width limits, keep pressure as the hydraulic screen pending a physical budget, and minimize normalized dissipation once flow feasibility is demonstrated. Use the verified adjoint for those responses. Further geometric expansion, area-rule adoption, and changes to the pressure budget or operating requirements remain approval decisions.

## Reproduction and figures

Sol job **2083835** is the reference run. The earlier **2083795** run produced the case, gradient, and figure files, but its launch wrapper passed a Windows carriage return to the final exit command and Slurm recorded a failure. The repeat used a Linux-format batch script and completed in 65 seconds. All case values and gradient histories matched the original run. Sol job **2083880** refreshed the figure labels using the same verified data. Both simulation runs remain in the project record. [Confirmed Slurm log](output/sol-2083835/slurm-2083835.out) · [Original Slurm log](output/sol-2083795/slurm-2083795.out) · [Run settings](output/sol-2083835/settings.json)

From the solver workspace on Sol, use the existing Brinkman Python environment. Set OpenBLAS and OpenMP to one thread each and request two CPUs so the two physical models can run in separate workers. The local archive path below supplies the same saved Re 50 probe used on Sol.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg
python -m project4_fem.transition_pilot --archive project4_fem/output/sol-2082990 --output pilot-output
python -m project4_fem.plot_transition_pilot --output pilot-output
```

| Figure | Contents | PNG | PDF |
| --- | --- | --- | --- |
| A | Original fillet region and expanded transition envelope | [PNG](output/sol-2083835/design-envelopes.png) | [PDF](output/sol-2083835/design-envelopes.pdf) |
| B | Measured side flow, target band, area failures, and pressure violations | [PNG](output/sol-2083835/measured-side-flow.png) | [PDF](output/sol-2083835/measured-side-flow.pdf) |
| C | Friction-only predictions and FEM width sweep | [PNG](output/sol-2083835/prediction-versus-fem.png) | [PDF](output/sol-2083835/prediction-versus-fem.pdf) |
| D | Reconstructed area versus measured side flow | [PNG](output/sol-2083835/area-side-flow.png) | [PDF](output/sol-2083835/area-side-flow.pdf) |
| E | Re 50 velocity magnitude with the reconstructed contour over the Brinkman solution | [PNG](output/sol-2083835/pilot-velocity-re-50.png) | [PDF](output/sol-2083835/pilot-velocity-re-50.pdf) |

[Pilot implementation](transition_pilot.py) · [Plotting implementation](plot_transition_pilot.py) · [Saved Re 50 FEM state](output/sol-2083835/velocity-re-50.npz) · [Master project handoff](../COOLANT_MANIFOLD_MASTER_HANDOFF.md)
