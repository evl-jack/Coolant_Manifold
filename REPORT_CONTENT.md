# Adjoint-based optimization of a coolant manifold junction

## 1. Summary
A sharp planar T-junction is reshaped using a 128-variable Brinkman density field and a custom Taylor–Hood finite-element solver at Re = 50. The main pipe and side outlet collar remain fixed. On the finest tested fitted-wall mesh (0.0625 mm target), the raw Navier–Stokes optimum reduces mechanical pumping power by **4.8372%**, from 9.275934 to 8.827238 µW/m depth, while increasing the side-flow fraction from 4.024206% to 8.464686%. The paired mesh screen passed. The serial adjoint sensitivity benchmark showed a 60.48× speedup over complete centered finite differences. The immersed endpoint passed its KKT audit. The nonlinear design also outperforms the Stokes-optimized geometry in the same nonlinear fitted-wall evaluator. Insert the live viewer at the top. Source: verification and nonlinear README/JSON; regenerate exact numbers using export_results.py.

## 2. Structure and formulation
A coolant manifold distributes pumped coolant. A sharp junction can alter local losses and branch delivery. This is a planar hydraulic study, not a thermal or validated 3D study. Coordinates in mm: main pipe [0,20]×[0,4]; side collar [9.5,10.5]×[6,7]; editable region [8,12]×[4,6]; original branch [9.5,10.5]×[4,7]. The fluid reference area is 83 mm². Density cells 0.25 mm, 128 independent unitless values m_i in [0,1], 1 fluid and 0 solid. The larger 960-entry array contains fixed values; do not call all 960 independent. The computational union is 89 mm². Fixed cells remain anchored through filtering/projection/reconstruction.

$$\min_m f=P(m)/P_{ref},\quad R(z,m)=0,\quad g=(\eta_{ref}-\eta(m))/0.001\le0,\quad 0\le m_i\le1.$$

The original immersed design replay at the same continuation stage sets both references; final beta=8 eta_ref=0.074147376839036. This is a team-selected nondegradation requirement, not an industry standard. The physical fitted reference flow is a separate validation metric, not the immersed constraint. No area or turning-angle constraint; those quantities are diagnostics. One independent PDE inequality is present; do not invent a second constraint.

$$P=-\int_{\Gamma_{ports}}(p+\tfrac12\rho|u|^2)(u\cdot n)\,ds,\quad \eta=Q_{side}/Q_{in},\quad R_{eff}=P/Q_{in}^2.$$

Here power is W per meter depth and flow is m²/s; R_eff is Pa·s/m². Keep outward normal and sign conventions consistent with source. Distinguish port power from dissipation and account for viscous boundary work and discrete divergence defect in the energy audit.

## 3. Analysis and optimization slide
Steady incompressible momentum: $\rho(u\cdot\nabla)u-\nabla\cdot(2\mu D(u)-pI)+\alpha u=0$, $\nabla\cdot u=0$. Parabolic prescribed inlet, no slip on walls, zero total traction at both outlets. Re=50; mean inlet speed 0.01254758565 m/s; inlet flow per depth 5.01903426e-5 m²/s. Confirm physical scales in code. P2 velocity/P1 pressure triangular elements, nonlinear Newton with the complete convection derivative $(\delta u\cdot\nabla)u+(u\cdot\nabla)\delta u$. Original analysis mesh 2,848 triangles / 13,359 DOFs. Use the actual weak residual, quadrature and elimination conventions in source.

The anchored Helmholtz mapping uses $H=M+r_f^2K$, $H_{EE}\tilde\gamma_E=M_E m_E-H_{EF}m_{0,F}$ for diagonal cell mass and fixed anchors; $r_f=0.5$ mm. Projection uses the archive's tanh Heaviside at threshold 0.5 (identity at beta=0). Brinkman interpolation $\alpha=\alpha_{max}q(1-\bar\gamma)/(q+\bar\gamma)$, alpha_max=100000 kg/(m³·s), q=0.1, beta=0,2,8. MMA iterates densities; physical contour at 0.5. Exact original raw field at start. The pilot used 100 unique-evaluation cap per stage; 97 total unique MMA evaluations reported. Read stage termination details instead of claiming every native stop converged.

**Slide diagram:** 128 densities → anchored filter → projection → Brinkman resistance → nonlinear FEM → port power & flow constraint → two discrete adjoints → MMA update → repeat. State variables u,p; objective and single inequality above; Re 50 is N2 nonlinear physics. Use real geometry and histories.

## 4. Verification
Use the existing ladder with exact expected values, computed values, error and pass/fail. Historical V1–V7 summaries are in historical draft, but verify provenance and applicability. Current archive contains mass/energy checks, analysis refinement, tangent/filter/transfer tests and gradient tests. Do not invent missing exact or independent hand-derived fixtures. Preserve original failure flags. See MISSING_DATA.md.

## 5. Sensitivities
From $R(z,m)=0$, $R_z dz/dm=-R_m$. For each scalar F, solve $R_z^T\lambda_F=\nabla_z F$ and compute $\nabla_mF=F_m-R_m^T\lambda_F$. Fixed-port direct term is zero. Reverse chain includes quadrature transfer, rational interpolation derivative $-\alpha_{max}q(q+1)/(q+\bar\gamma)^2$, projection derivative, and anchored filter transpose $H_{EE}^T\psi=s$, $\nabla_mF=M_E^T\psi$. Use ADJOINT_METHOD.md for full derivation of both responses. Forward + two adjoint RHS for both gradients.

15/15 directional FD tests passed. Original Taylor order screen passed 9/15, with six flagged; larger-step follow-up passed all six (orders 2.0010–2.0060). The initial JSON's aggregate passed=false is intentional and must not be overwritten. Complete 128-variable centered gradient (257 forward solves including baseline) infinity errors: 2.376e-12 normalized power and 1.099e-12 side fraction. Serial same-node benchmark 120.379 s FD vs 1.990 s adjoint, speedup 60.48×. Scaling fixed at 52,347 analysis DOFs: 32/128/512 controls, 10.672/12.472/12.741 s adjoint and 40.339/112.083/395.261 s FD elapsed with eight workers. Distinguish aggregate work from parallel elapsed and same-node serial comparison.

## 6. Optimization and KKT
Final normalized objective 0.992919627; g=-9.97890521; side-flow slack 0.997891 pp. Projected stationarity 1.369e-8 vs tolerance 1e-4; primal, dual, complementarity residuals zero. 49 lower, 73 upper, 6 interior density variables. Flow multiplier 0; inactive constraint implies zero first-order local cost for a sufficiently small increase in required flow. For box bounds use $w=\nabla f+\lambda\nabla g$, $\|m-\operatorname{clip}(m-w,0,1)\|_\infty$. KKT is numerical first-order stationarity only, not global optimality or fitted-wall certification. Show the actual archived history and residual map.

## 7. Nonlinearity: Stokes versus Navier–Stokes
Removing convection from momentum and its Jacobian makes the fixed-density Stokes state linear; keep the same inlet, design mapping, common final side floor and objective including kinetic transport. Both endpoints pass KKT. Stokes satisfies the immersed constraint under its own evaluator, but fails under nonlinear evaluation (6.990948% vs 7.414738% floor). Navier–Stokes optimum gives 8.412628% under nonlinear immersed evaluation. Under identical nonlinear fitted-wall evaluation, Stokes delivers 1.948290% side flow and 9.105665 µW/m, while NS gives 8.464686% and 8.827238 µW/m. The NS optimum consumes 3.06% less power than Stokes and increases side share by 6.516397 pp. Paired descriptive uncertainties 0.013075 µW/m and 0.022710 pp; all three fitted geometries pass mesh screens. The Stokes geometry fails to preserve the original fitted side flow. The raw density mean absolute difference is 0.61866, symmetric-difference area 2.078946 mm². Stokes area 86.228145 mm², turn 88.247971° vs NS 85.506975 mm², 54.678079°. The observed total times 220.520 s Stokes/217.466 s NS reflect different evaluation counts (300/97) and environments, not intrinsic solver speed. Embed real nonlinear figures and complete six-row cross table from JSON.

## 8. Results and numerical qualification
| Geometry | Power µW/m | Resistance Pa·s/m² | Side % | Area mm² | Wall turn ° |
|---|---:|---:|---:|---:|---:|
| Original fitted | 9.275934 | 3682.284 | 4.024206 | 83.000 | 90.000 |
| Raw NS fitted | 8.827238 | 3504.165 | 8.464686 | 85.507 | 54.678 |
| Smoothed fitted (provisional) | 8.858017 | 3516.383 | 8.229312 | 85.458 | 55.558 |
| Stokes fitted (NS evaluator) | 9.105665 | 3614.692 | 1.948290 | 86.228 | 88.248 |

All rows finest 0.0625 mm; use regenerated JSON for display. Raw NS physical saving 0.448696 µW/m vs paired uncertainty screen 0.011197 µW/m, combined energy residual 0.001788 µW/m. Mesh criteria two successive changes <=0.5% power and <=0.020 pp side flow. Smoothed failed the coarse-medium side change 0.02505 pp; retain as provisional. Immersed primary paired screen also passes but finite resistance is not an impermeable wall. The original immersed flow differs substantially from original fitted flow. The nominal Brinkman penetration length ~0.10010 mm has ~1.60 finest target intervals; avoid claims of resolved zero leakage. No matched-area/matched-flow control available. Increased area and changed outlet allocation mean smoothing alone cannot be credited for the entire saving. No imposed 45° limit in this model and measured maximum exceeds 45°.

## 9. Bonus, 10. limitations and AI disclosure
N2 nonlinear Newton + implicit differentiation is a candidate bonus, subject to actual nonlinear benchmark evidence. Density topology with anchored filter/projection is another candidate, but required comparable fixed-topology control remains missing unless independently documented. No global optimum, no validated 3D, no heat-transfer or turbulence claim. No rigorous numerical error bound. AI disclosure must be completed by team: actual agent-written code, human-derived checks and a documented bug caught through verification. Do not invent student contributions. See historical draft for one exploratory SST issue, not automatically a bug in the graded laminar solver.
