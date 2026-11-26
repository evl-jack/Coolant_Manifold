# Coolant manifold mathematical optimization demonstration

**Decision.** The verified discrete adjoint predicted a local reduction in energy, and independent forward solves confirmed the prediction. The derivative-driven SLSQP attempt reduced the excessive side flow but did not produce a feasible or lower-energy constrained design. The fitted-wall study distinguishes the optimizer endpoint from the separate adjoint-descent diagnostic.

## Original manifold and fixed optimization problem

Use the original 20-control mirrored fillet with \(a_i\in[0,1]\) and \(a_0=(0.5,\ldots,0.5)\). The controls displace the wall on \(y=5.0\)–\(6.5\) mm about the centerline \(x=10\) mm. The planar inlet at \(x=-1\) mm and main outlet at \(x=21\) mm are each 2.00 mm tall. The side outlet at \(y=11\) mm is 1.25 mm wide. The original map keeps the throat at 1.00 mm and fixes straight guides and ports. The existing P2–P1 incompressible laminar FEM uses Re 50, a unit-flow parabolic inlet, zero natural outlet traction, and no slip on the fixed walls. The fixed optimization mesh has **2,116 triangles**. The later 26-control pilot and metering reach have different geometry and area scopes. This mathematical demonstration uses the original approved objective, design family and screens. [Original formulation](FORMULATION.md) · [Scope](OPTIMIZATION_DEMONSTRATION_SCOPE.md)

Let \(S(u)=\frac12(\nabla u+\nabla u^\mathsf T)\) and let \(R(u,p,a)=0\) be the assembled momentum and continuity system with Brinkman drag \(\alpha(a)u\), where \(\alpha=20{,}000(1-\gamma_h(a))^2\). The objective is \(J=(D_\mathrm{viscous}+D_\mathrm{Brinkman})/111.1367471732773\), where the two terms integrate \(2S:S\) and \(\alpha|u|^2\) in the nondimensional solver. Side flow is \(f_s=Q_s/Q_\mathrm{in}\). Static pressure loss is inlet area-mean pressure minus the signed-flow-weighted outlet pressure. The area variable is diffuse solid area \(A_s=\int_\mathrm{strips}(1-\gamma_h)\,dA\) over the complete original strips.

| Quantity | Optimization rule |
| --- | --- |
| Normalized dissipation | Minimize \(J\) |
| Side flow | \(0.195\le f_s\le0.205\), with both inequalities scaled by 0.005 |
| Static pressure loss | At most **95.31983921402909**, scaled by that cap |
| Diffuse solid area | At most **0.396620614844318 mm²**, scaled by that cap |
| Original shape bounds and clear width | \(0\le a_i\le1\), width \(\ge1.00\) mm by construction |

These are enforced screens for this numerical problem. The pressure cap comes from 105% of the original Brinkman reference, and a hardware-derived allowance remains outstanding. The course permits a laminar inertial comparison for a manifold. Its final optimized inertial versus non-inertial comparison remains open. [Archived limits](output/sol-2082990/limits.json) · [Official Project 4 specification](https://designinformaticslab.github.io/DesignOptimization2025/project4.html)

## Discrete adjoint and independently measured predictions

With \(w=(u,p)\), solve \(R_w^\mathsf{T}\lambda=J_w^\mathsf T\), then form \(\nabla_a J=J_a-R_a^\mathsf T\lambda\). Geometry enters the design-mesh residual through drag, with \(\alpha_{a_i}=-40{,}000(1-\gamma_h)\gamma_{h,a_i}\). The cost also has a direct drag derivative. Side fraction and static pressure use the exact discrete flux and trace functionals, and diffuse area has a direct density derivative. The constraint row scalings and signs are applied to the reduced Jacobian. The established adjoint reuses the factorized state Jacobian for the different responses. [Solver and adjoint](fillet.py) · [Registered optimizer protocol](optimization_demo/PROTOCOL.md)

The independent Sol gate **2085839** compared objective, flow, pressure and area responses in three directions using 42 forward states and seven centered step sizes from 0.04 to 0.000625. All **12** derivative and Taylor comparisons passed. The largest finest-step error was **0.432% of its permitted allowance**. The historical tolerance \(\max(10^{-7},0.002|D|)\), observed Taylor-order criterion, \(10^{-11}\) free residual and \(10^{-9}\) relative mass gate were retained. This evidence covers the immersed 20-control formulation. [Raw verification and endpoint states](optimization_demo/sol-2085871/gradient/gradient-gate.json) · [Independent protocol](optimization_verification/PROTOCOL.md)

At \(a_0\), the independently evaluated test direction was \(-\nabla J/\|\nabla J\|_\infty\). Each forecast was saved before a cold-start PDE solve. Negative change is less energy. The flow changes in the last column are measured percentage-point increases.

| Maximum coefficient change | Predicted \(\Delta J\), % of original | Measured \(\Delta J\), % of original | Measured minus predicted | Pressure change, predicted / measured | Side-flow change (pp) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.0400 | −0.049542 | −0.049392 | +0.000150 | −0.068800 / −0.068376 | +0.028926 |
| 0.0100 | −0.012385 | −0.012377 | +0.000009 | −0.017200 / −0.017179 | +0.007344 |
| 0.0025 | −0.00309635 | −0.00309603 | +0.00000033 | −0.00430000 / −0.00429920 | +0.001841 |

All three perturbed PDEs passed the state gates. For the 0.04 step the predicted and measured diffuse solid-area changes were **−0.0093840** and **−0.0093853 mm²**. The energy prediction error relative to the forecast fell from **0.3040%** to **0.0105%** as the step shrank. Each forward solve still missed the flow band. [Prediction figure](optimization_demo/sol-2085871/math/predicted-measured.png) · [All predicted and measured constraints](optimization_demo/sol-2085871/math/prediction-table.csv)

## Genuine constrained SLSQP attempt

Sol job **2085871** started the original problem at \(a_0\) with exact adjoints, all four original nonlinear inequalities, the box bounds, and no relaxation. SLSQP used its internal line search with the established 35-iteration maximum and `ftol=1e-9`. It evaluated **68 distinct converged PDE states**, used **18** gradient-base points and recorded **50** other trial points. SciPy returned its iteration counter at **22** and status **8**, a positive directional derivative in the line search. There were no PDE gate failures. Major callbacks, bases and other trial evaluations remain separately labeled. [Objective versus accepted iteration](optimization_demo/diagnostics/accepted-objective.png) · [Trace figure including callbacks](optimization_demo/sol-2085871/math/objective-history.png) · [Every evaluation and constraint](optimization_demo/sol-2085871/math/iteration-table.csv) · [Full prospective line-search and derivative trace](optimization_demo/sol-2085871/math/events.jsonl) · [SLSQP termination](optimization_demo/sol-2085871/math/termination.json)

| Immersed L1 result | Original | SLSQP returned point |
| --- | ---: | ---: |
| Dissipation, viscous plus Brinkman | 111.13675 | 112.26549 |
| Objective relative to original | 1.000000 | **1.010156** |
| Side / main share of inlet | 25.1011% / 74.8989% | **24.2187% / 75.7813%** |
| Static pressure loss | 90.78080 | 92.33097 |
| Diffuse solid area (mm²) | 0.39662061 | 0.39881273 |
| Minimum width (mm) | 1.00000 | 1.00000 |

The returned state left **3.7187 percentage points** of excess side flow, exceeded the diffuse area cap by **0.0021921 mm²**, and retained **2.9889** pressure units below the cap. Dissipation rose by **1.12875**, or **1.0156%**, while side share fell **0.8824 points**. Seventeen controls were within the existing \(10^{-5}\) activity tolerance of a box bound. The terminal objective-gradient norm was **0.01145**. The endpoint remained outside the feasible set, so stationarity and KKT multipliers for a feasible optimum are unavailable. Slurm job completion confirms the script exited and its results were saved. The reported optimizer status governs the mathematical outcome. [Returned design and raw responses](optimization_demo/sol-2085871/math/candidate.json)

The terminal gradient base has a local constraint conflict. Using its archived Jacobian and the true coefficient bounds, a linear program minimized the linearized side share while satisfying the linearized lower-flow, pressure and area inequalities inside a coefficient step radius of 0.04. The best linearized side share was **24.2215%**, slightly above the measured endpoint. At radius 0.0025 even linearized area restoration was unavailable. This identifies a geometry and area tradeoff in the terminal linearization, consistent with the outer line-search failure and converged inner solves. The diagnostic supplies neither an infeasibility proof for the nonlinear local neighborhood nor a global bound. [Mathematical results and terminal diagnostic](optimization_demo/RESULTS.md) · [Independent numerical replay](optimization_verification/VERIFY.md)

## Independent fitted-wall PDE verification

Sol job **2085900** rebuilt the original, endpoint and prespecified 0.04 descent contours and independently solved the PDEs with no-slip fitted walls. Each geometry used **2,116**, **8,464** and **19,044** triangles and separate Stokes → Re 25 → Re 50 startup. All nine Re 50 states and all continuation states passed residual and conservation checks. Maximum residual was **8.14×10⁻¹²** and maximum relative mass imbalance **1.47×10⁻¹³**. Pressure, side flow and dissipation passed the established L2-to-L3 checks for every geometry. The descent comparison was triggered because the hard SLSQP run returned a larger objective. The fitted integral measures viscous work with the Brinkman term absent. [Fitted solver and exact mesh inputs](optimization_verification/PROTOCOL.md) · [Raw fitted comparisons](optimization_verification/sol-2085900/fitted/fitted-check.json)

| Finest fitted Re 50 result | Original | SLSQP endpoint | 0.04 adjoint descent |
| --- | ---: | ---: | ---: |
| Side share of inlet | 24.88809% | **23.93718%** | 24.91770% |
| Main share of inlet | 75.11191% | 76.06282% | 75.08230% |
| Side-flow excess above 20.5% (pp) | 4.38809 | **3.43718** | 4.41770 |
| Static pressure loss | 91.12968 | 92.78176 | **91.05744** |
| Pressure headroom | +4.19016 | +2.53808 | +4.26240 |
| Viscous dissipation | 111.40037 | 112.59288 | **111.34198** |
| Dissipation change | — | **+1.19251 (+1.07047%)** | **−0.05839 (−0.05241%)** |
| Analytic geometric fluid area (mm²) | 51.903208 | 51.901597 | 51.912695 |
| Analytic area change (mm²) | — | **−0.001611** | **+0.009486** |
| Finest polygonal fluid area (mm²) | 51.903556 | 51.901933 | 51.913042 |
| Minimum passage width (mm) | 1.000 | 1.000 | 1.000 |

The optimizer endpoint reduced fitted side share by **0.95091 points** and increased fitted energy and pressure. It also lost physical fluid area. The fitted gain of the descent diagnostic moved in the predicted energy direction and was positive on all three meshes. From L2 to L3 its paired energy saving changed by **0.000171** versus a fine saving of **0.058390**. Side flow increased by **0.02962 points**. This independently confirms that the adjoint direction produced a lower-dissipation contour, with a flow-split penalty. It is a separate measured result from the hard-constrained attempt. All three fitted contours were pressure- and width-compliant at Re 50 and above the side-flow band. [Reconstructed contours and fitted results](optimization_verification/sol-2085900/fitted/fitted-geometry-responses.pdf) · [Three-mesh results table](optimization_verification/sol-2085900/fitted/fitted-comparison.csv)

With convection disabled, the same fine original and endpoint contours gave **27.4210%** and **26.8373%** side share, a **0.5837-point** paired reduction. Both Stokes pressures exceeded the common cap. These are fixed-shape comparisons of the influence of inertia. The project still lacks two feasible optima for the course's optimized nonlinear comparison.

## Cross-code context and research decisions

The previously verified Re 50 OpenFOAM comparison corroborated the flow mechanism on related fitted geometry. It measured a paired side-flow reduction of **1.1467 points**, compared with **1.1765** in fitted FEM. The separately sampled metering shapes showed a design-dependent pressure discrepancy. OpenFOAM's fixed outlet pressure and zero normal velocity-gradient condition differs from FEM's natural zero total traction. The earlier absolute dissipation discrepancies are about **0.25%–0.30%** and do not bound the error in the new paired energy change. A pressure-only reclassification cannot move the hard-optimizer endpoint into the target flow band. The new fitted states make a new targeted OpenFOAM run unwarranted for the registered decision. The small energy gain is supported by an independent fitted-wall re-solve and remains untested in a second code or hardware. [Archived OpenFOAM refinement](next_openfoam/OPENFOAM_VALIDATION.md) · [Independent CFD review](optimization_cfd/REVIEW.md)

The mathematical connection is verified locally for the approved immersed geometry. The hard-constraint run reveals a tradeoff between side-flow reduction, the area bound and the energy objective. A broader geometry or changed hardware requirement would need approval and gradient verification for the actual new parameterization. This finite run leaves feasible designs elsewhere open. The measured energy result and feasibility must remain separate in the final presentation.

## Reproduction and independent oversight

| Sol phase | Job | Workers | Elapsed | Allocated CPU-hours | Measured CPU-hours | Peak Slurm batch RSS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Independent derivative gate | 2085839 | 4 | 11 s | 0.01222 | 0.00776 | 424,548 KiB |
| Prediction and hard SLSQP | 2085871 | 1 | 33 s | 0.00917 | 0.00824 | 178,792 KiB |
| Fitted mesh comparison | 2085900 | 2 | 90 s | 0.05000 | 0.04435 | 1,972,616 KiB |

Each worker used one BLAS and OpenMP thread. SLSQP itself was serial. Total cost was **0.07139 allocated CPU-hours** and **0.06034 measured CPU-hours**. The base Git revision was **917c4600f7d95b28a19257110eea59e16d55e1ca**. The stage inputs are identifiable by the archived source and per-file digests. The fitted job used a supplemental source archive and the exact output vectors written by the mathematical run. [Math sources](optimization_demo/sol-2085871/source.tar.gz) · [Math and gradient hashes](optimization_demo/sol-2085871/input-sha256.txt) · [Fitted input manifest](optimization_verification/sol-2085900/inputs.sha256) · [Independent Sol accounting and reproduction](optimization_hpc/SOL_REPRODUCIBILITY.md)

Independent compliance reviewed the frozen objective and limits, the derivative gate, prospective prediction logging, acceptance classification, endpoint selection and conditional fitted diagnostic. Optimization, numerical verification and compliance used separate source ownership. HPC and CFD reviewed on available slots. Optimization and verification were launched with GPT‑6.2 Astra and ultra requested, compliance and CFD with GPT‑6.2 Astra and high requested. HPC inherited its settings. Exact lead model and effort were unavailable to this session. [Independent compliance record](optimization_compliance/REVIEW.md)
