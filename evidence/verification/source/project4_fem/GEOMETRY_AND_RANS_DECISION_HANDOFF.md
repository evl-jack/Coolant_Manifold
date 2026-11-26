# Project 4 decision handoff — geometry feasibility and RANS adjoint

**Decision.** The current 26-control geometry has a measurable effect on flow split. At Re 50, the finest fitted-wall mesh gives **24.8871%** side flow at the baseline and **23.7106%** for the screened candidate. The **1.1765-percentage-point** reduction exceeds the **0.1228-point** empirical uncertainty in the paired gain. The candidate still misses the upper flow target by **3.2106 points**. The next geometry worth reviewing is a smooth metering reach in the side branch, paired with explicit area compensation. Keep Re 50 as the laminar reference and plan the RANS adjoint for a separate turbulent operating case.

## Pilot audit

Job **2083835** remains the pilot record. Both rank checks returned 26, all 26 sampled states converged, and both physics models passed the four directional adjoint tests using the original centered-difference and Taylor criteria. The two density-map derivative errors were below **3.8×10⁻⁹**. At the finest derivative step, the worst error consumed **0.57%** of the allowed tolerance in Stokes and **1.50%** at Re 50. These results concern the laminar residual and its responses. [Pilot handoff](TRANSITION_PILOT_HANDOFF.md) · [Case table](output/sol-2083835/pilot-cases.csv) · [Stokes derivatives](output/sol-2083835/gradients-re-0.json) · [Re 50 derivatives](output/sol-2083835/gradients-re-50.json)

| Item | How it was handled | Result or limitation |
| --- | --- | --- |
| Geometry | Enforced by the map and verified on the FEM mesh | Bounded controls, mirrored displacements, continuous connected passage, fixed interfaces and joints, and minimum width of 1.00 mm. |
| Operating conditions | Prescribed in the state solve | Unit parabolic inlet, no-slip fixed walls, zero natural outlet traction, and common evaluation planes. |
| Area | Computed from the analytic contour and density quadrature, then screened against the matched baseline | Nine shapes passed both measures per model. A scalar geometric compensation restored area in the combined probes. This rule was approved for the pilot. |
| Flow and pressure | Numerically screened against 19.5%–20.5% and 95.3198 | No shape passed the flow band. All Re 50 samples passed pressure and all Stokes samples exceeded the cap. |
| Numerical checks | Residual and mass gates, plus gradient tests | Pilot residual ≤9.73×10⁻¹⁴ and relative mass imbalance ≤1.38×10⁻¹³. These gates certify the discrete solution and baseline derivatives. |
| Other responses | Reported as diagnostics | Dissipation, side/main ratio, and nominal-solid leakage. The pilot had no leakage limit or optimizer stationarity claim. |

The historical 20-control Brinkman baselines were **27.6483%** in Stokes and **25.1011%** at Re 50. The 26-control pilot baselines were **27.9587%** and **25.3151%**, shifts of **+0.3104** and **+0.2140 points**. The expanded mesh and the new diffuse-wall interval both affect the comparison. The matched gains below use the same mesh and wall treatment within each pair. The sampled cases establish local evidence and leave feasibility across the full design space open. [Historical recovery audit](RECOVERY_DECISION_HANDOFF.md) · [Historical formulation](FORMULATION.md)

## Matched refinement and wall-treatment results

Sol job **2084088** completed the preregistered 24 solves in 76 seconds using four workers. Both shapes used the same saved 26 coefficients, with mesh levels of **2,500**, **10,000**, and **22,500** triangles per method. The fitted meshes use no-slip walls along the piecewise-linear approximation of the analytical density midpoint. The Brinkman meshes use the same subdivisions and the archived penalty 20,000. The port planes at x = −1 mm, x = 21 mm, and y = 11 mm and the opening widths are unchanged. All states passed the port, nonlinear residual and mass checks. Maximum nonlinear residual was **4.10×10⁻¹³** and maximum relative mass imbalance was **7.99×10⁻¹³**. [Prerecorded protocol](output/sol-2084088/MATCHED_VERIFICATION_PROTOCOL.md) · [All matched results](output/sol-2084088/matched-verification.csv) · [Sol log](output/sol-2084088/slurm-2084088.out)

| Model and wall treatment, finest mesh | Side flow, base → candidate | Gain (pp) | Static pressure loss, base → candidate | Total dissipation, base → candidate | Solid RMS speed, base → candidate |
| --- | ---: | ---: | ---: | ---: | ---: |
| Stokes, Brinkman | 27.9481% → 26.9560% | 0.9921 | 107.9158 → 108.8642 | 107.2554 → 108.2048 | 0.7212% → 0.8712% |
| Stokes, fitted | 27.4197% → 26.2853% | 1.1344 | 108.4889 → 109.5566 | 107.8292 → 108.8976 | Outside the fluid domain |
| Re 50, Brinkman | 25.3072% → 24.2502% | 1.0570 | 90.4486 → 91.9802 | 110.9124 → 111.9458 | 0.9068% → 1.0350% |
| Re 50, fitted | 24.8871% → 23.7106% | 1.1765 | 91.1308 → 92.8386 | 111.4011 → 112.5361 | Outside the fluid domain |

Pressure and dissipation are nondimensional. Pressure loss is area-mean static inlet pressure minus signed-flow-weighted static outlet pressure. Total dissipation includes viscous work and, for Brinkman, penalty work. On the fine Re 50 candidate those contributions were **111.5440** and **0.4018**, respectively. Solid RMS speed uses points with γ below 0.05 and is expressed relative to mean inlet speed. Both shapes retain 1.00 mm width. The candidate adds **0.000820 mm²** of reconstructed fluid area. The finest Brinkman and fitted meshes integrate **+0.000802** and **+0.000816 mm²** relative to their own baselines. All area screens passed.

On every paired level 2-to-3 comparison, changes stayed within **0.0077 points** in side flow, **0.059%** in pressure, **0.026%** in dissipation and **0.88%** in nominal-solid velocity where defined. The wall treatments missed the declared absolute flow and pressure agreement criteria for both shapes and both physics models. At Re 50 the fitted wall lowered side flow by **0.4201 points** at baseline and **0.5397 points** for the candidate and raised candidate pressure loss by **0.925%**. Increased mesh resolution did not remove that wall-treatment discrepancy at this penalty. The gain envelopes were **0.1453 points** in Stokes and **0.1228 points** at Re 50, each smaller than either method’s measured gain. The geometry effect is resolved for the tested pair. All Stokes pressure values exceed the cap, all Re 50 values pass it, and every matched case misses the flow band. [Full differences and screening flags](output/sol-2084088/matched-verification.json) · [Refinement figure PDF](output/sol-2084088/matched-side-flow.pdf)

![Measured side-flow fractions across the three matched FEM mesh levels](output/sol-2084088/matched-side-flow.png)

## Geometry with enough potential authority to review

The tested area-neutral curvature mode changed Re 50 side flow by **−0.0335 points**. Full coordinated transition narrowing alone changed it by **−0.7921 points**, or **−0.6234 points** when compensated from the baseline fillet. The short editable transition therefore limits the tested resistance change. The fitted candidate leaves **3.2106 points** to reach the band. Using the original pilot’s **3.7258-point** gap as a second calibration gives the same ordering of alternatives. [Measured pilot cases](output/sol-2083835/measured-side-flow.pdf)

The screening calculation uses slit resistance proportional to $R_b=\int_{5}^{11} w(y)^{-3}\,dy$, holds effective main-path resistance fixed, and estimates the new split as $f_s=r/(r+R_b^\mathrm{new}/R_b^\mathrm{seed})$ with $r=f_s^\mathrm{seed}/(1-f_s^\mathrm{seed})$. The fitted seed needs **20.5%** more branch resistance to reach the upper band and **24.3%** to hit 20%. The original pilot seed needs **24.0%** and **27.9%**, respectively. These are friction screens. Junction losses, outlet traction effects and changes in effective driving pressure require FEM. The prior Re 50 friction estimate overstated the effect of full transition narrowing by **0.3633 points**. [Friction proxy compared with FEM](output/sol-2083835/prediction-versus-fem.pdf)

| Candidate region | Shape and compensation | Analytical result and tradeoff |
| --- | --- | --- |
| **Preferred for review: lengthen the narrow reach downstream** | Release just y = 8–10.5 mm inside the current guide. Use a 1.00 mm width from y = 8.5–10.0 mm, with a 0.50 mm cubic ramp on each end and tangent joins. Keep the last 0.50 mm at 1.25 mm plus both main interfaces. Symmetric cuts suit a constant-depth milled or extruded passage, subject to cutter radius and stock clearance. | **+27.4%** friction resistance before compensation, predicting **19.62%** from the fitted seed and **20.07%** from the pilot seed. Loses **0.5000 mm²**. Even fully opening the current fillet leaves **0.2523 mm²** of compensation to find. Enlarging a low-resistance branch shoulder and tuning the plateau length are candidates for review. Full compensation will change the estimate. |
| **Lower-scope alternative: use existing compensation only** | Open the current fillet to its allowable maximum and shorten the narrow reach to about **0.49 mm** between the same ramps. Keeps the minimum width and needs no new shoulder widening. | Area-neutral in the analytical geometry. **+10.6%** resistance gives **21.93%** from the fitted seed or **22.42%** from the pilot seed. This particular compensated combination falls short in the proxy. |
| **Keep the guide: extend the neck into the upper fillet** | Smoothly narrow the upper fillet toward 1.00 mm over y = 5.6–6.5 mm and compensate in the broader shoulders. Maintain tangency and check curvature for the actual cutter and wall stock. | Requires about **0.0803 mm²** of gross area redistribution. **+2.4%** compensated resistance gives **23.29%** from the fitted seed or **23.80%** from the pilot seed. Curvature remains a secondary tuning option. |

In the friction screen the long reach can span the target, but the area deficit and pressure headroom determine whether it is worth building. The current fitted candidate has **2.4812** pressure units of headroom. A fixed-main linear resistance estimate would need **96.75** pressure units merely to reach 20.5%, above **95.3198**. That estimate scales the seed loss by main flow share and cannot establish violation of the FEM static pressure metric. Review whether there is room for a smooth compensation pocket at the branch shoulder first. If approved, begin with a few fitted-wall candidates that preserve the 1.00 mm width, port geometry and target, and screen pressure and area before any optimization.

Baseline-relative area is compatible with preserving total coolant inventory at constant depth and a fixed exterior envelope. A minimum area bound does allow inventory to increase and says nothing about wall stock, local throat tolerance, bend radius, tool access, or allowable coolant volume. Obtain the manufacturing interpretation before adopting a production rule. Continue to treat 19.5%–20.5% as a hard hydraulic constraint and 95.3198 as the provisional pressure cap from the historical Re 50 baseline. Its hardware basis is open. Geometric feasibility, hydraulic feasibility, primal convergence, KKT stationarity and physical validity each require their own evidence. KKT evaluation is deferred until there is a specified, feasible optimization problem and a computed candidate optimum.

## RANS model decision

The course lists Stokes or Navier–Stokes with Brinkman as a manifold project example, with inertial nonlinearity around Re 10–200. Its RANS guidance requires differentiation of the turbulence model wherever RANS is used. The Re 50 laminar reference is physically appropriate for the present reduced case, and its verified laminar adjoint remains applicable to that case. No evidence supports sustained turbulence at this operating point. For perspective, the ideal wide slit would give a hydraulic-diameter Reynolds number near 100 when based on \(D_h\approx 2H\). Pipe-transition data provide context for the scale, with a different geometry and inlet disturbance level. Confirm whether the course team expects a separate RANS deliverable for this project. [Project 4 requirements](https://designinformaticslab.github.io/DesignOptimization2025/project4.html) · [Avila et al. pipe-transition experiment](https://doi.org/10.1126/science.1203223)

| Repository component | Available today |
| --- | --- |
| [Project 4 Taylor–Hood FEM](fem.py) and [fillet adjoint](fillet.py) | P2 velocity, P1 pressure, convection, Brinkman drag and transpose solve for drag-dependent design gradients. There are no turbulence state variables or SST residual blocks. |
| [Original 2D finite-volume solver](../rans_solver/solve.py) and [wall reconstruction](../rans_solver/brinkman.py) | A segregated SIMPLE/SST-style forward implementation, wall distance and optional near-wall immersed treatment. Useful as a cross-check after matching closures and boundary data. It has no verified coupled FEM adjoint. |
| [Custom 3D solver](../rans3d/breakdown.md) and [Sol OpenFOAM case](../sol_rans_3d/README.md) | Separate SST-style and SST forward workflows. They supply a higher-speed water operating precedent and external comparisons. The meshes, pressure conventions and validation status differ from this FEM study. |

Plan a **separate water operating point** in a demonstrably turbulent regime, using measured flow, actual depth and hydraulic diameter, inlet development and roughness. As a planning example only, 2 m/s water in a 2 mm wide-slot gap would give \(Re_H\approx 4{,}000\) and \(Re_{D_h}\approx 8{,}000\) in the wide-slit limit. Establish turbulent inlet conditions using measured profiles or a validated developed-flow precursor. Preserve the normalized mass flow for paired designs at the new operating point and specify a pressure budget for that hardware regime separately. Reconsider the planar model if depth or secondary flow changes the split.

Use **wall-resolved Menter SST-2003m** as the first closure to verify. Its strain-based eddy-viscosity and production limiters address adverse pressure gradients and separated shear layers, with two additional scalar fields. Start on fitted no-slip walls and keep the exact model variant and constants consistent across references. Use the inlet \(k\) and \(\omega\) profiles from measurements or a precursor. For bounded sensitivity trials, intensity \(I\) and length scale \(\ell\) can estimate \(k=1.5(UI)^2\) and \(\omega=\sqrt{k}/(C_\mu^{1/4}\ell)\), with \(C_\mu=0.09\). Where \(k=0\) at walls, prescribe the variant’s finite near-wall \(\omega\) from first-point spacing and molecular viscosity. Target first off-wall \(y^+\le1\), resolving the separation zone and repeating at finer spacing. Use zero diffusive flux on outflow for both scalars and specified positive turbulence data on any inflow portion of an outlet. Monitor backflow and validate the treatment at the fixed pressure planes. [SST equations and boundary prescriptions](https://tmbwg.github.io/turbmodels/sst.html) · [Wall-spacing study](https://tmbwg.github.io/turbmodels/flatplate_val_ypluseffect.html)

## Fully coupled adjoint implementation and verification plan

Use \(w=(u,p,k,\omega)\) with residual \(R=(R_m,R_c,R_k,R_\omega)\). Add continuous P1 scalar spaces to the existing P2–P1 velocity-pressure spaces. Check scalar convection stabilization, positive scalar states and stiff near-wall source terms on fitted meshes before introducing design-dependent diffuse walls. Start with a segregated or continuation-based initialization and drive the same coupled discrete residual to convergence with damped Newton. Compute the final Jacobian from that residual. The partial-derivative structure on a fixed incompressible mesh without pressure-dependent scalar boundary data is:

$$
R_w =
\begin{bmatrix}
M_u & G & M_k & M_\omega \\
D & 0 & 0 & 0 \\
K_u & 0 & K_k & K_\omega \\
W_u & 0 & W_k & W_\omega
\end{bmatrix}.
$$

- \(M_u\) includes mean convection, strain dependence of eddy viscosity, and drag on immersed meshes. \(G\) and \(D\) are the discrete pressure and divergence couplings. \(M_k,M_\omega\) include eddy stress, the active SST limiter branches and pressure/traction dependence when required by the convention.
- \(K_u,W_u\) include advection, signed boundary convection, strain production and stabilization. \(K_k,K_\omega,W_k,W_\omega\) include turbulent diffusivity, blending, cross-diffusion, production, destruction and positivity regularization. In SST-2003m the isotropic stress can enter modified pressure. When used, recover static pressure \(p=p^*-2\rho k/3\) in dimensional variables, differentiate this recovery, and apply a consistent outlet stress convention.
- Assemble analytic velocity-pressure terms already verified. Use element-local float64 automatic differentiation of the same scalar, stress, stabilization and boundary kernels used in the primal. Use a differentiable fitted mesh map for \(R_d\), including quadrature transforms, wall distance, wall spacing and boundary data. For a later Brinkman version, also differentiate drag, turbulence damping and the reconstructed wall distance with respect to \(\gamma(d)\). Verify that the chosen distance model and solid treatment converge toward the fitted result. Record limiter switches. Where smoothing is needed, verify the smoothed model’s bias and differentiate exactly that model.
- Solve \(R_w^\mathsf{T}\lambda=J_w^\mathsf{T}\) and compute \(dJ/dd=J_d-\lambda^\mathsf{T}R_d\). Use separate right-hand sides for side fraction, static pressure loss and the approved dissipation response. Specify the RANS objective before optimization. Mean-flow energy loss from molecular plus eddy stress includes transfer to turbulent kinetic energy. Its budget must account for that transfer and transported \(k\) without counting the same energy loss twice. Frozen turbulence can serve as a sensitivity diagnostic only. [Discrete adjoint for turbulent topology optimization](https://orbit.dtu.dk/en/publications/topology-optimization-of-turbulent-flows/)

Proposed gates for the future RANS verification, declared here before implementation:

| Check | Acceptance and evidence to retain |
| --- | --- |
| Coupled primal | Scaled free-DOF infinity norm ≤10⁻⁹ separately for momentum, continuity, k and ω, with fixed documented scales. Relative mass imbalance ≤10⁻⁸. Verify turbulence equations’ integrated flux/source balances to ≤10⁻⁶ of their absolute term sums. Verify wall and inlet data, scalar signs and outlet backflow handling. |
| Jacobian | Centered directional differences of each scaled residual block for isolated and mixed state directions, and of \(R_d\) for geometry directions. Require relative defect ≤10⁻⁶ with a 10⁻⁸ absolute floor over a stable step range, plus a normalized transpose-product defect ≤10⁻¹⁰. Record active limiters and repeat near switches with the chosen smooth treatment. |
| Adjoint | Re-converge both perturbed states in two random directions, a collective direction and a transition-only direction. Use the seven pilot steps, 0.04 through 0.000625, and compare each normalized response separately. Require final error ≤max(10⁻⁷, 0.002 \|directional derivative\|) and Taylor order between 1.7 and 2.3 over a usable interval or final remainder <10⁻¹². |
| Resolution and validation | Across three meshes, require medium-to-fine side share within 0.10 points, pressure and dissipation within 1%, and tested directional gradients within 5% with a 10⁻⁶ floor in normalized units. Verify wall shear and profiles on a turbulent channel or flat plate against an independent implementation of the same variant. Validate separating-flow reattachment, velocity and stress against the backward-facing-step data, then compare manifold split and pressure against hardware or a matched 3D reference. Report SST versus Wilcox k–ω closure sensitivity with compatible boundary inputs and compare the spread against hydraulic margins. |

The fine laminar meshes have **102,933 state degrees of freedom**. Two P1 scalars add roughly **23,000** degrees of freedom on the same connectivity, before wall refinement. Budget **1–6 CPU minutes per converged RANS state** as a planning range, then benchmark assembly, factorization, memory and iterations on Sol. Four directions and seven centered steps cost up to **56 primal solves**, roughly **1–6 CPU hours per design per grid** at that rate. Run independent perturbations in parallel. For a converged state the cost of each adjoint is a transpose solve reusing one numerical factorization or a compatible preconditioner. Separate effort is needed for wall-distance and boundary validation. [Channel benchmark](https://tmbwg.github.io/turbmodels/channelflow_val.html) · [Backward-facing-step measurements](https://tmbwg.github.io/turbmodels/backstep_val.html)

## Decisions and presentation

1. Approve whether to release y = 8–10.5 mm of the guide for a metering land and identify the allowable shoulder area, exterior stock, cutter radius, and dimensional tolerance. Retain the terminal opening and bore interfaces.
2. Confirm the final area rule and obtain the hardware pressure allowance at the declared port planes. Decide whether the Stokes comparison is only a reduced model for the course or must meet the same hydraulic screens.
3. Confirm whether the course expects a separate RANS deliverable. Approve a turbulent operating point, matching inlet data and validation reference before implementing the fully coupled model. Work stops at this decision report.

| Slide | Use the verified evidence |
| --- | --- |
| 3 | End the workflow at screened geometry and matched verification. Keep constrained optimization as a future stage. |
| 4 | Label Stokes and laminar Navier–Stokes at Re 50 as the two solved equations. Show SST as a separate planned regime. |
| 5 | Show \(\alpha=20{,}000(1-\gamma)^2\), the analytic contour, and what nominal-solid leakage measures. |
| 6 | Show the coupled \((u,p,k,\omega)\) residual and transpose solve. Mark RANS and its adjoint as proposed for Project 4. |
| 7 | Show the 20+6 control envelope, fixed ports and three matched meshes. Label the area formulation as a pilot screen. |
| 8 | Replace equal-outlet uniformity with the 19.5%–20.5% side share and the defined static pressure metric. |
| 9 | State the fixed flow, two zero-traction outlets, separate wall treatments, and mesh and wall sensitivity. |
| 10 | Correct the presentation configuration to **one inlet and two outlets**. Use the verified 2D T-geometry and converged Re 50 field. |
| 11 | Use the matched refinement plot, paired gain and uncertainty, and the verified 26-control derivative table. Distinguish convergence from target satisfaction. |

The presentation mapping uses the workspace’s [Slide 3–11 outline](../../presentaion.md). Additional figures are available for the [26-control geometry](output/sol-2083835/design-envelopes.pdf), [pilot area tradeoff](output/sol-2083835/area-side-flow.pdf), and [converged Re 50 velocity field](output/sol-2083835/pilot-velocity-re-50.pdf). The matched figure was generated on Sol from job 2084088’s state records.

Reproduce with the existing Sol Python environment. Set `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1` in a four-CPU allocation.

```bash
python -m project4_fem.matched_verification --pilot project4_fem/output/sol-2083835 --output matched --workers 4
python -m project4_fem.plot_matched_verification matched
```
