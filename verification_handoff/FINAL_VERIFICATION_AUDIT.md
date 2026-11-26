# Final verification audit — frozen 128-variable fixed-pipe pilot

Evidence base: bba3589e359daa6489b13446155b677bffd31b91. Official [rubric](https://designinformaticslab.github.io/DesignOptimization2025/project4.html#rubric). This is a source/archive audit, not a new numerical campaign. No website, geometry, solver or optimization changed. Classification: final-model = exact pilot family or its archived independent replay; historical = other mathematical verification; exploratory = different geometry/model/settings; failed = prescribed gate missed; incomplete = unavailable/inconclusive evidence. See [raw tables](VERIFICATION_TABLES.json), [index](EVIDENCE_INDEX.md), [limitations](REMAINING_LIMITATIONS.md).

## 1. Final formulation
Planar per unit depth, water rho=998.2 kg/m³, mu=.001002 Pa s, H=.004 m, U=Re mu/(rho H), Re=50, Q'=UH=5.01903426e-5 m²/s. Main rectangle 20×4 mm, collar x=9.5–10.5, y=6–7 mm; editable rectangle x=8–12, y=4–6 mm, 128 independent .25 mm cells. Original sharp T area 83 mm²; analysis envelope 89 mm². Fixed no-slip walls, inlet u=6y(1-y) nondimensional, zero symmetric total traction at outlets. Viscous pressure scale mu U/H, drag scale H²/mu, power scale mu U² (W/m). P2/P1 Taylor–Hood, 25-point mapped Gauss quadrature, sparse Newton with residual-decreasing line search. 
Nondimensional equations:
```math
Re(u\cdot\nabla)u-\nabla\cdot(2S(u))+\nabla p+\alpha^*u=0,\quad\nabla\cdot u=0.
```
S=(grad u+grad u^T)/2. Weak form integral [2S(u):S(v)+Re v·(u·grad)u+alpha* u·v-p div v-q div u], less loads/traction. Tangent has both convection terms Re[(delta u·grad)u+(u·grad)delta u], plus viscous, pressure, continuity and drag blocks. Fixed DOFs eliminated; 13,359 total, 12,571 free, 788 prescribed. [FEM](../evidence/verification/source/project4_fem/fem.py), [model](../evidence/verification/source/project4_fem/project43_topology/fixed_pipe_pilot_model.py).
Actual NLP min f=P/P_ref subject to g=(eta_ref-eta)/.001<=0, 0<=m_i<=1. References are original immersed state at each stage; beta8 P_ref=8.893614 µW/m and eta_ref=.074147376839. No independent area or turning cap, no second PDE-performance constraint in final pilot. Disclose potential formulation-rubric gap; the floor is baseline-derived, not a hardware standard. Port objective P=-mu U² sum over ports integral [(p+Re|u|²/2)(u·n)] ds with outward normals; eta=Q_side/Q_in with inlet positive. R_h=P/(Q')². Do not substitute static pressure loss or total traction work.

## 2. Historical FEM verification ladder
Raw [baseline verification](../historical/baseline_verification.json) and [verification implementation](../evidence/verification/source/project4_fem/verify.py). All nine compound archive flags pass; split here into tests. Nondimensional results, L2 errors as implemented. Historical execution job/hash not supplied.
| Test | Analytical expectation | Actual / measured error | Prescribed gate | Status |
|---|---|---|---|---|
| Interpolation | P2 nodal delta | 0 | 1e-12 | Pass |
| Partition of unity/gradient | sum N=1, sum grad N=0 | both 0 | 1e-12 | Pass |
| Quadrature | integral x^i y^j=i!j!/(i+j+2)! | max moment error 1.110223e-16 | 1e-12 | Pass |
| Poiseuille | Q=1 | Q=.9999999999999918, relative error 8.216e-15; velocity/pressure L2 1.184e-14/7.054e-14 | 1e-8 | Pass |
| Rotated channel | rotation .61 rad, Q=1 | Q=.9999999999999974; relative error 2.554e-15; velocity/pressure L2 9.313e-15/3.917e-14 | 1e-8 | Pass |
| Global mass | Qout=Qin | max 3.477e-14 across straight/rotated/uniform Brinkman | 1e-10 | Pass |
| MMS velocity | P2 smooth expected order 3 | n=4/8/16 errors .0525773/.00601501/.000702917; orders 3.12780/3.09714 | last >2.6 | Pass |
| MMS pressure | P1 smooth expected order 2 | errors .0378547/.00404450/.000377283; orders 3.22644/3.42224 | last >1.7 | Pass |
| Nonlinear MMS | analytic Re=35 field | n=12 velocity L2 .001701212; pressure .001025404 | velocity <.01 | Pass |
| Newton tangent | centered Jd at h=1e-5 | rel 1.076192e-11 | 1e-7 | Pass |
| Uniform Brinkman | Q=.290570054755 | .290569429301, rel 2.152507e-6 | 5e-4 | Pass |
| Layered alpha100 | total .409744853735 | .409744859224, rel 1.33938e-8 | 5e-4 | Pass |
| Layered alpha400 | total .241589204240 | .241589243195, rel 1.61244e-7 | 5e-4 | Pass |
| Layered alpha1600 | total .173374318848 | .173374495586, rel 1.01940e-6 | 5e-4 | Pass |
| Layered solid flux | .08653698293/.02679814315/.00712465950 | .08653714370/.02679905513/.00712705078; rel 1.86e-6/3.40e-5/3.36e-4; fractions .211198/.110928/.041108 decrease | 1e-3 and monotone | Pass |
| Benchmark adjoint | 24 Gaussian controls, 2 responses, 2 directions | smallest-step max rel 1.32962e-7; min Taylor slope 1.99898 | 2e-4; slope >1.8 | Pass |

Mathematical purpose: basis tests catch interpolation/assembly and integration mistakes. Poiseuille u=6y(1-y), p=12(2-x) in 2×1 nondimensional channel, exact outlet shear traction; rotation tests gradient and normal transformations. Mass tests flux integration and continuity. MMS streamfunction psi=8f(x)f(y), f(t)=t²(1-t)², u=(psi_y,-psi_x), p=x+y/2, forcing=-laplacian u+grad p+Re(u·grad)u. Exact prescribed trace and pressure pinned at origin. Nonlinear residual sequence .0163031592→.00235551742→1.43605e-7→2.68584e-16; supports fast local Newton decrease, not a fitted asymptotic rate at roundoff. Uniform Brinkman -u''+alpha u=G, G=12, alpha=25, u=0 at walls; u=G/alpha[1-cosh(sqrt(alpha)(y-.5))/cosh(sqrt(alpha)/2)], Q=G/alpha[1-2tanh(sqrt(alpha)/2)/sqrt(alpha)]. Layered core |y-.5|<=.25 parabolic, outer strips hyperbolic, u,u' continuous, exact flux integrated; these test finite leakage, not impermeability. Benchmark adjoint alpha=25(1-gamma)^2, pressure-difference and upper-half outlet flux; h=.04,.02,.01,.005,.0025. Historical data includes exact gradients and FD values, not final mapping.

Source identity: repository tree shows identical fem.py Git blob 2826e76dc590d8f2d6d23c09f52a6f684cdffc4f across pilot, verification and nonlinear source copies; verify.py blob 69bf2d05f553774d593bbf0951cc9781a977c2a2. Historical report has no original executed source hash or job ID, so historical executable identity not established. No verified established Sol runtime in this session; no replay. Targeted replay command (not executed): from archived source parent, python -m project4_fem.verify --output NEW_DIRECTORY/verification.json, unchanged tolerances.

## 3. Complete final density chain and implicit derivatives
Editable E, fixed F. Raw expansion fixes binary anchors. Anchored Helmholtz A f_E=M_E m_E-H_EF m_F plus eliminated top-halo forcing; A includes halo diagonal, H=M+r_f²K, r_f=.5 mm. f_F=m_F. Forward A^-1 M_E d, reverse M_E^T A^-T s. Anchor exchange means no global density conservation claim. Projection threshold t=.5, gamma=[tanh(beta t)+tanh(beta(f-t))]/[tanh(beta t)+tanh(beta(1-t))], derivative beta sech²(beta(f-t))/denominator, beta0 identity. Rational drag alpha=alpha_max q(1-gamma)/(q+gamma); derivative -alpha_max q(q+1)/(q+gamma)^2; q=.1, alpha_max=1e5 kg/(m³ s). Cell constant resistance at each assigned quadrature point, alpha*=H² alpha/mu. Transfer transpose sums quadrature cotangents. Source assertions preserve raw/filtered/projected fixed values and zero fixed-fluid drag; cell area check rtol 1e-12, fixed area 81 mm², widths 4/4/1 mm. Reconstruction threshold .5 clipped to design rectangle, no change outside and exact fixed pipe/collar coverage asserted. Analytic rank of linear editable map is 128 (A nonsingular, M positive diagonal); positive projection and nonzero drag derivative retain rank at finite beta. This is a mathematical inference, not numerical shape-map rank; thresholded contour is nonsmooth. Earlier measured rank 20 belongs to another map. See [fixed model](../evidence/verification/source/project4_fem/project43_topology/fixed_pipe_pilot_model.py), [anchor parent](../evidence/verification/source/project4_fem/project43_topology/wall_displacement_density_v1.py), [adjoint method](../evidence/verification/source/project4_fem/project43_topology/verification/ADJOINT_METHOD.md).

For free state z, R(z,m)=0 implies dz/dm=-R_z^-1 R_m. For F=f or eta, R_z^T lambda=F_z^T and total derivative F_m-R_m^T lambda. Fixed ports/prescriptions imply F_m=0, lambda extended by zero on fixed DOFs. Quadrature drag cotangent -w_q u_q·lambda_q; reverse transfer and H²/mu, rational derivative, projection, anchored transpose and compact extraction. If boundary values were design-dependent this would need a lifting derivative; they are fixed here.
For h=p+Re|u|²/2, delta P=-mu U² sum integral [(delta p+Re u·delta u)(u·n)+h(delta u·n)] ds; f_z=P_z/P_ref. eta_z=(Qside_z-eta Qin_z)/Qin, g_z=-eta_z/.001. P2 cubic kinetic flux degree six; four-point line quadrature exact on straight edges; independent five-point port check.
Final original, gray, beta0/2/8, random, junction, collar directions: interior centered; bound inward second-order (-3F0+4Fh-F2h)/(2h). Gates transfer 1e-12, port parity 1e-10, state 1e-10, mass 1e-9, residual derivative 1e-6, scaled adjoint 1e-4, Taylor order >1.8. Independently decomposed gray and beta8 state tangents complete error <=7.22e-16 and centered relative <=4.81e-11; prescribed values preserved. Final isolated filter/projection/transfer/rank metrics not all independently recovered; end-to-end checks and source assertions support chain, historical component tests distinct.
15/15 FD directions pass. Initial Taylor 9/15; six flags beta_0/junction_cell, beta_0/collar_cell, beta_2/junction_cell, beta_2/collar_cell, beta_8/junction_cell, beta_8/collar_cell. Insufficient consecutive above-floor intervals at original steps (floor 5e-11); supplement uses .04/.02 and original steps, no gate relaxation; all six pass, reported orders 2.0010–2.0060. Preserve both plots and full records. Full coordinate centered gray-field 128 test: max absolute errors 2.376e-12 (normalized power), 1.099e-12 (side fraction), 257 forward solves; not centered coordinate testing at bound endpoint. Source [original checks](../evidence/verification/adjoint-evidence/adjoint_merged.json), [supplement](../evidence/verification/adjoint-coordinate-supplement.json).

## 4. Adjoint cost (Sol 2408534)
| Case | Controls | DOFs | Adjoint elapsed s | FD elapsed s | Aggregate FD CPU s | FD forward solves incl baseline |
|---|---:|---:|---:|---:|---:|---:|
| pilot-128 |128|13359|1.990|120.379|119.795|257|
| scale-32 |32|52347|10.672|40.339|202.702|65|
| scale-128 |128|52347|12.472|112.083|755.298|257|
| scale-512 |512|52347|12.741|395.261|2999.272|1025|
Serial 60.48× uses unrounded times (rounded ratio 60.49). Scaling uses eight FD workers on identical FEM mesh, not pilot mesh; elapsed ratios not serial speedups. CPU work distinct from elapsed and aggregate wall work. One forward state, one factorization, two transpose solves; each FD state supplies both responses. Setup/import excluded. [Timing method](../evidence/verification/source/project4_fem/project43_topology/verification/TIMING_METHOD.md).

## 5. MMA/KKT
Binary sharp start, continuation q=.1 alpha=1e5 beta0→2→8, 100 unique eval cap each, previous selected warm start; no restoration/restart in final pilot.
| beta | unique | native termination | normalized f | side % | area mm² |
|---|---:|---|---:|---:|---:|
|0|17|design tolerance|.992965509|7.977077|87.149976|
|2|18|design tolerance|.993157106|8.084346|86.592271|
|8|62|objective tolerance|.992919627|8.412628|85.506975|
97 unique, 99 objective callbacks; local summary elapsed 223.773 s (stage sum 217.466); not Sol. Stage-specific normalizations not one fixed objective. Independent KKT job 2408529 replayed archive. L=f+lambda g, lambda>=0, primal g<=0 and bounds, complementarity lambda g=0, stationarity ||m-clip(m-grad L,0,1)||inf. Lower bound grad L>=0, upper <=0, interior zero. Final g=-9.97890521, slack .997891 pp, lambda=0, projected stationarity 1.369e-8, primal/dual/complementarity zero. Lower/upper/interior 49/73/6. Gates 1e-6 primal, 1e-8 dual, 1e-4 stationarity, 1e-5 complementarity. Zero multiplier means zero local marginal cost of small floor tightening while slack remains. First-order local stationarity, not global optimum/second-order proof. [KKT](../evidence/verification/kkt-results/fixed_pipe_kkt.json), [pilot](../evidence/verification/input/summary.json).

## 6. Physical evidence
| Case | h mm | triangles | P µW/m | resistance Pa s/m² | side % |
|---|---:|---:|---:|---:|---:|
|immersed original|.25|2848|8.893614|3530.514|7.414738|
|immersed original|.125|11392|8.892091|3529.910|7.431694|
|immersed original|.0625|45568|8.892020|3529.881|7.437211|
|immersed raw|.25|2848|8.830644|3505.517|8.412628|
|immersed raw|.125|11392|8.831105|3505.700|8.424696|
|immersed raw|.0625|45568|8.831774|3505.965|8.429885|
|fitted original|.25|4132|9.282894|3685.047|4.038784|
|fitted original|.125|16466|9.277632|3682.959|4.027304|
|fitted original|.0625|65970|9.275934|3682.284|4.024206|
|fitted raw|.25|4338|8.834575|3507.077|8.443024|
|fitted raw|.125|16986|8.828640|3504.721|8.459634|
|fitted raw|.0625|67860|8.827238|3504.165|8.464686|
|fitted smoothed|.25|4236|8.867440|3520.124|8.196856|
|fitted smoothed|.125|16961|8.859554|3516.993|8.221901|
|fitted smoothed|.0625|67865|8.858017|3516.383|8.229312|
Raw fitted saving .448696 µW/m (4.8372%), +4.440481 pp side; descriptive paired uncertainty .011197 µW/m, combined absolute raw closure .001788 µW/m. Both successive mesh changes pass .5% power/.02 pp side. Smoothed saving 4.5054%, but first split change .02505 pp fails .020; provisional. Original/raw/smoothed areas 83/85.507/85.458 mm², max turn 90/54.678/55.558° measured at .25 mm secant. Raw recommended. Area and flow allocation change; no causal isolation of rounding. 
Immersed fine reduction .060246 µW/m (.6775%), +.992674 pp, uncertainty .001984 µW/m, raw closure sum .002994 µW/m; optimizer-grid .708037%. Finite resistance, gray filter, threshold geometry and no-slip representation differ. Fine immersed original side 7.437211% vs fitted 4.024206%; raw 8.429885 vs fitted 8.464686. Fine immersed original gross ingress/egress ~4.316/4.287% of inlet, net .0299%; net near zero does not establish impermeability. Penetration .10010 mm, finest h gives 1.60 nominal intervals; not penalty convergence. Energy: report advective port power, viscous work, dissipation, drag dissipation, weighted divergence, raw/corrected closure and discrete weak reaction closure separately. Global mass conservation not pointwise divergence-free. Individual state/mass/energy/leakage metrics and job IDs in [mesh raw](../evidence/verification/mesh-evidence/mesh_report.json) and machine-readable handoff.

## 7. Matched Stokes/NS nonlinear comparison
Remove convection from residual/tangent, retain Re=50 inlet, same map, port energy, common side floor 7.414738%. See [nonlinear report](../evidence/nonlinear/nonlinear_component.md).
| Design | evaluator | P µW/m | side % | meets floor |
|---|---|---:|---:|---|
|original|Stokes|10.200410|7.545020|yes|
|original|NS|8.893614|7.414738|yes|
|Stokes optimum|Stokes|10.117120|7.414742|yes|
|Stokes optimum|NS|8.907379|6.990948|no|
|NS optimum|Stokes|10.279525|9.204333|yes|
|NS optimum|NS|8.830644|8.412628|yes|
Stokes fitted NS replay h .25/.125/.0625: power 9.113994/9.106854/9.105665, side 1.947624/1.941523/1.948290%. Finest resistance 3614.6924 vs NS raw 3504.1646. Stokes-minus-NS .278427 µW/m, -6.516397 pp; uncertainty .013075 µW/m and .022710 pp; combined raw closure .002065 µW/m. Not equal allocation. Stokes area 86.228145 mm², turn 88.247971°, mean abs density difference .61866, symmetric difference area 2.078946 mm². Stokes KKT stationarity 1.069e-5, complementarity 4.789e-8, multiplier .0010054; 300 unique evaluations vs NS 97, timing not pure per-solve speedup. Inertia changes preferred shape and feasibility; no turbulent or 3D validation.

## 8. Historical campaign ledger
- Pre-FEM: supplied Project43_searchable.txt, COOLANT_MANIFOLD_PRE_FEM_HANDOFF.md. Cartesian nonuniform mesh connectivity, 20×4 mm 80-cell rectangle; SIMPLE/Rhie–Chow/upwind/SST, 2 m/s, rectangle ~102 and T ~420 iterations. Exploratory density drag (1-gamma)^2, budget changes and 3-control pinched geometry: fitted .5/.25 mm side 28.06/28.60%, pressure 1784/1520 Pa. Not final.
- Older 3D OpenFOAM mesh failures 2053128/2057077, repaired 2059820/2059943 pass mesh only; custom square duct error .819→.206%. URANS 2061743 17.42% side, 7.83 kPa window 10–20 ms; pressure drift 7.67→7.99 kPa. Exploratory. Symmetry sweep 2079036 29/59 converged, max-side 2079305 31/65 and 15.63% side. Fine 2081608 momentum defect 39–43× gate despite stable hydraulics; failed. Resolved sensitivities only 2079711/2079893.
- Axisymmetric annular manifold, SU2 and earlier paired leakage investigations: different cylindrical integration and nW, historical leads; exact raw sources/job records not independently recovered here. Do not transfer to final planar W/m.
- [20-control mirrored fillet](../evidence/verification/source/project4_fem/MATHEMATICAL_OPTIMIZATION_HANDOFF.md): 2085839 all 12 derivative/Taylor comparisons pass, 2085871 SLSQP status8, 68 converged states, infeasible and not lower energy; independent descent h=.04 predicted -.049542%, measured -.049392%, misses split band. Historical gradient support, failed optimizer. [26-control transition](../evidence/verification/source/project4_fem/TRANSITION_PILOT_HANDOFF.md) different geometry/area, exact endpoint not recovered.
- Planar baffle/shoulder 2349590 targets 50%; uploaded report.json at alpha1e6, .125mm: immersed balanced/revised side 46.3100/47.7182%, power 56.0865/35.7600 µW/m, gross ingress 21.114/16.469%. Exploratory, different constraints.
- [Gate2](../evidence/verification/source/project4_fem/project43_topology/gate2_evidence.md): 2395055 serialization, 2395294 corrected runner, 2395378 finer one-sided sweep passed. Earlier Neumann filter transpose 1.78e-15, conservation 1.22e-16, rank20 on 240/960/3840 cells, derivative 5.12e-9 vs1e-6; transfer 2.28e-16 vs1e-12, reduced/full 3.40e-16. Not anchored map. 2395413 both SLSQP/MMA 250-budget endpoints fail KKT (stationarity diagnostic1); feasible trials not converged optimization. 2395642 native infeasible-start analytical MMA violation .940979 despite callback FD 8.43e-11; feasible start recovered (.5,.5) within1.97e-9. 2395874 Phase I/II recovered within1.74e-9 and 18 KKT fixtures pass. Lesson: restore feasibility, independent endpoint checks.
- [Gate3](../evidence/verification/source/project4_fem/project43_topology/gate3_evidence.md): fitted fixtures 2395603 pass, binary penalty screen 2395875 all 9 combos fail at least one gate per fixture. At alpha1e7 sharp split .15227 pp vs .02, obstructed .11730 pp and power error1.08358% vs .5. Resolved local layer 2396361 still split .168779/.129576 pp. Failed physical parity, not necessarily coding error.
- SST: uploaded Project43_SST_Newton_searchable.txt: 2353228 bounded iteration stagnation sweep217; 2354828 transpose <=2.61e-15 but failed directional gate 1e-4 at two amplitudes; mixed k 2.13e-4 at1e-7; branch-aware stencil far point unchecked. Native prescribed-state 2377694/2377706 controller-reported leads, not retrieved independently. No converged turbulent design or new runs.
- Earlier 1mm displacement 2406318 beta4 cap500, beta8 not started; reported fitted diagnostic 9.36281%,9.31143 µW/m vs immersed10%,8.78525, turn66.44° vs45°; incomplete, not final.
- Final fixed-pipe local pilot + independent Sol replays 2408522–2408661: final-model evidence.

## 9. Rubric mapping
| Category | Points | Evidence/qualification |
|---|---:|---|
|Formulation|10|§1 actual objective/bounds, one PDE inequality only, engineering target rationale|
|Analysis and verification|20|§2 historical exact/MMS/balance/invariance/nonlinear ladder; source identity gap|
|Sensitivity|20|§3–4 both responses, full chain, flags and follow-up, cost|
|Optimization|15|§5 history and independent KKT|
|Nonlinear|15|§7 Stokes vs NS and benchmark|
|Website/3D|10|website agent eight requirements; visualization not CFD|
|Reproducibility|5|raw archive/scripts/jobs; clean-checkout single command untested|
|Clarity|5|handoff + limitations; student AI statement confirmation|
|N2 bonus|+5|strong evidence, not promised points|
|Topology bonus|+5|partial, no compatible optimized fixed-topology comparator found|
|Other bonuses|—|no qualifying evidence established|
