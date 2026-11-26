# Gate 2 verification on ASU Sol

Sol scratch reports are under `/scratch/jtfoste6/p43-topology-study-20261125/`. The analysis is planar at Re = 50 with a fixed sharp T-junction baseline. The penalty of 100,000 kg/(m³·s) in these verification fixtures is provisional.

| Job | Scope | Outcome |
|---|---|---|
| 2395055 | First common verification | Density and three baselines completed. Serialization of a NumPy boolean interrupted the report runner. |
| 2395196 | Isolated dependency staging replay | Installed the pinned NLopt, geometry, and test packages in a job-local path. This run used the original report runner. |
| 2395294 | Corrected report runner | Density, baseline, reconstruction and SLSQP smoke checks completed. All interior derivative directions met tolerance. The one-sided sweeps ended too early to get two passing residual differences per direction. The pinned MMA binding lacked an optional interface used by the wrapper. |
| 2395378 | Finer one-sided sweep | All discrete-flow verification checks passed using unchanged thresholds. |
| 2395413 | Matched optimizer replay | Both algorithms reached the shared 250-evaluation budget. Each sampled a feasible design. Neither passed the independent endpoint KKT audit. |

## Density and state verification

| Design cells | Wendland basis rank | Condition number | Smallest and largest singular values |
|---:|---:|---:|---:|
| 240 | 20 | 2.20813 | 1.21304, 2.67854 |
| 960 | 20 | 2.25943 | 2.36858, 5.35164 |
| 3,840 | 20 | 2.26810 | 4.71575, 10.69581 |

The exact binary initial geometry and 83 mm² area matched at all three spacings. The filter's largest transpose error was 1.78e-15, largest conservation error 1.22e-16, and observed refinement orders 1.986 and 1.996. The largest final-step density derivative error was 5.12e-9 against a 1e-6 limit. All singular values, coverage and finite-difference tables are in the Sol JSON reports.

On the 808-triangle, 3,909-DOF flow verification mesh, the transfer dot-product error was 2.28e-16 against a 1e-12 limit, and three nontrivial reduced-to-full chain comparisons agreed within 3.40e-16. Mass imbalance was at most 1.12e-15 and nonlinear residual at most 9.40e-12, against 1e-9 and 1e-10 limits. The minimum reported interior Taylor order was 1.955. At the binary bounds, extending the sweep recovered 3 passing steps for full topology and 4 for the reduced model. The smallest-step residual errors were 1.72e-8 and 7.90e-9 against a 1e-6 limit. Corresponding largest response or constraint directional errors were 1.21e-6 and 1.72e-6 against a 1e-4 limit.

## Fixed-fluid seams and baseline comparison

The analysis spacing was 0.25 mm and the design spacing was 0.50 mm. At the two side-seam design cells, projected density was 0.631662 and resistance was 574.04 kg/(m³·s). The fixed collar retained unit fluid density and zero resistance. The three baseline comparisons showed material sensitivity to both finite resistance and filtering.

| Baseline | Triangles | Side flow | Port power (µW/m) |
|---|---:|---:|---:|
| Sharp fitted | 2,656 | 4.03949% | 9.28043 |
| Binary immersed, unfiltered | 3,232 | 5.86163% | 8.86624 |
| Binary-initialized immersed, filtered | 3,232 | 8.57790% | 8.14703 |

Python 3.12.16, NumPy 2.5.4, SciPy 1.18.1, NLopt 2.9.1, Shapely 2.1.2 and Triangle 20250106 are recorded in job 2395294. The ten contour and constrained-meshing reconstruction fixtures passed. The SLSQP analytical smoke test returned (0.5, 0.5) to 1.4e-15 and passed the independent KKT audit. The 80-evaluation SLSQP pilot found a feasible trial but hit its budget and returned a non-stationary, infeasible last major iterate. The matched optimizer replay is the follow-up to that attempt.

## Matched optimizer results

Both algorithms used an exact binary start with 8.548799% side flow and a common power reference. The best feasible SLSQP trial was 4.729675 P/P0. At the evaluation cap the last major iterate was infeasible, with maximum scaled g of 424.1070. SLSQP used 251 forward solves and 502 adjoint right-hand sides including the terminal audit. The MMA/CCSA endpoint was also its best feasible trial at 4.726301 P/P0 and a maximum scaled g of −4.1153e-5. It used 250 forward solves and 500 adjoint right-hand sides. The independent LP found no valid multiplier witness for either point, with diagnostic projected stationarity 1 in both cases. MMA returned design-tolerance termination on the analytical constrained problem from an infeasible initial point while violating one constraint by 0.940979. That failure requires isolation before optimization verification can close.

Sol job 2395642 isolated the analytical MMA failure. Independent callback finite differences agreed within 8.43e-11 versus the 1e-8 comparison tolerance. Native scalar constraints, native vector constraints and the project adapter returned the same infeasible zero-start result. From the feasible start (0.4, 0.4), all three returned (0.5, 0.5) within 1.97e-9 and passed the independent KKT audit. This identifies the known-problem failure as native infeasible-start behavior of the pinned algorithm. The production handling of infeasible initial designs remains pending review.

## Two-phase recovery on Sol 2395874

The [reviewed restoration protocol](restoration_protocol.md) recovered the analytical optimum from the original infeasible origin. Phase I reached (0.325, 0.325) in three unique evaluations and both Phase II optimizers reached (0.5, 0.5) within 1.74e-9 with zero projected stationarity residual. All 18 independent KKT cases passed on Sol.

The full-topology 240-cell Phase I restored the original binary design to 50.095% side flow and 50.9850 mm² using 12 forward solves, 24 adjoint right-hand sides and 34 Newton steps in 3.72 seconds. The separate audit of that common start took one forward solve and two adjoint right-hand sides. The power at the feasible start was 6.047511 P/P0.

| Phase II optimizer | Forward solves | Adjoint RHS | Newton steps | Wall time (s) | Best feasible P/P0 | Endpoint |
|---|---:|---:|---:|---:|---:|---|
| SLSQP | 250 | 500 | 751 | 81.61 | 4.729430 | Evaluation limit. Last major iterate violated the split constraint by 408.209 scaled units. |
| MMA/CCSA | 250 | 500 | 610 | 73.54 | 4.730613 | Evaluation limit. Feasible endpoint and no KKT witness. |

The separate last-major-iterate audit added one forward solve and two adjoint right-hand sides for SLSQP. The MMA audit reused the cached terminal evaluation. Both audits found a diagnostic projected stationarity residual of 1 after the multiplier LP failed to find a witness. Neither endpoint is verified stationary.

The 20-variable reduced Phase I reached only 19.283955% side flow in its best observed trial, with approximately 62.4304 mm² fluid area. The 0.05-margin attempt consumed nine forward solves and the zero-margin fallback consumed 32. This is a local result at the provisional verification fixture. The run supplies no global infeasibility proof and no feasible reduced Phase II start. Phase I costs and inequalities are archived. Trial power logging was added to source after this run and is absent from its Phase I history.

## Reduced design-space diagnostic on Sol array 2396219

The reduced diagnostic maximized side flow subject to area and box constraints in the same provisional verification setting. Eight deterministic starts all gave 19.283955% side flow at the same box corner within roundoff, with 62.43037 mm² fluid area. Each passed the independent KKT audit for this diagnostic objective. The largest projected residual was 8.89e-16. The eight solves used 221 forward solves in their optimizer loops and provided no evidence that the current basis can reach 49.9% side flow at this fixture. The result does not establish global infeasibility.

A bounded least-squares fit of the full-space feasible anchor gave raw-density RMSE 0.33268 and projected-density RMSE 0.26412. A separate forward replay gave 18.66185% side flow in the fitted reduced density, compared with 50.09500% for the original full-space anchor. The exact box-only pointwise bound arrays and density coverage maps are archived beside the numerical reports. They precede application of the area constraint.

## Full-topology restart diagnostic on Sol array 2396219

Both optimizers restarted from their best feasible sampled checkpoint from job 2395874. Each restart used a fresh optimizer instance and the original constraints. A separate forward replay and directly integrated port responses supported the endpoint KKT audit.

| Optimizer | Total Phase II evaluations | Final checkpoint side flow | Area (mm²) | Final checkpoint P/P0 | Projected stationarity | KKT |
|---|---:|---:|---:|---:|---:|---|
| MMA/CCSA | 500 | 49.900007% | 70.830526 | 4.729402 | 6.14688e-5 | Pass |
| SLSQP | 1,000 | 49.901312% | 68.927384 | 4.728689 | Approximately 1 | Fail |

The independently replayed MMA power was 38.49628 µW/m. Its primal and dual residuals were zero. The LP bracket for best admissible projected stationarity was [6.14613e-5, 6.14688e-5], and complementarity was 1e-5 within the prescribed criterion. The matched SLSQP endpoint was feasible and exhausted its diagnostic budget. Its independently bracketed minimum projected stationarity was [0.99999905, 1].

Both inward directional checks passed. The centered checks at the archived endpoints gave power and split Taylor slopes above 1.94 and directional errors within their thresholds. The aggregate flagged the absence of a Taylor slope for the affine β=0 area response. Its first-order Taylor remainder was already at roundoff. Source now checks the affine remainder against the established roundoff floor and retains the separate finite-difference error requirement. Sol 2396915 replays that correction at both archived endpoints.

The centered MMA residual check in 2396915 had one passing perturbation before truncation dominated the next larger step. The follow-up retained the original sequence and added finer β=0 perturbations. Sol 2397105 passed both inward and centered probes at both checkpoints with unchanged thresholds. The smallest centered residual-action errors were 2.61e-9 for MMA and 8.56e-10 for SLSQP. The separate stationarity results in the table are unchanged.
