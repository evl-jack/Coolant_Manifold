# Numerical verification status

**Rungs 1–4 pass. Rungs 5 and 6 fail on incomplete verification.** Both of the final rungs produced results consistent with the analytical reference for the checks performed. Each still needs one part of its verification demonstrated numerically.

This evaluation compares the [analytical reference](verification_results.md) with the [numerical results](Numerical%20verification.md). The reports did not set acceptance tolerances before the runs. The limits below are the criteria used for this retrospective review. The verdicts cover the cases run by the verification solver in `fem.py`.

| Rung | Verdict | Evidence and criterion |
| --- | --- | --- |
| 1. Exact strip solution | **PASS** | Flow and driving power agree with the closed-form solutions within 0.000007%, and peak velocity within 0.000018%. The review limit for these quantities is 0.01%. Wall shear is 0.159% low in the hollow case and 0.197% low with Brinkman drag active, within the 0.25% review limit. |
| 2. Mesh convergence | **PASS** | The finest refinement step gives orders of 2.979 for velocity, 2.000 for pressure, and 1.992 for the velocity gradient. Each is within 0.1 of its expected order of 3, 2, and 2, respectively. |
| 3. Global balance | **PASS** | The outward end fluxes cancel to numerical precision. The maximum force imbalance is 0.0254% of the driving force, below the 0.1% review limit. Relative power imbalance is about $10^{-13}$, below the $10^{-8}$ review limit. |
| 4. Rotation | **PASS** | After mapping the rotated solution back, the largest nodal differences are $9.0\times10^{-19}\ \mathrm{m/s}$ in velocity and $6.0\times10^{-18}\ \mathrm{Pa}$ in pressure. Both are below $10^{-12}$ in their respective units. |
| 5. Nonlinear benchmark | **FAIL — incomplete** | The manufactured fields and Jacobian action agree with their references. Newton converged in two full steps, reaching the numerical floor before the run could establish a local quadratic convergence rate. |
| 6. Response gradients | **FAIL — incomplete** | The fixed-driving-force flow and power tangents agree with the hand sensitivities and centered differences. The depth Taylor remainder shows second-order scaling. The numerical report contains no adjoint comparison. |

## Rung 5. Measure the local Newton convergence rate

Newton reduced the residual from $3.17\times10^{-6}$ to $2.46\times10^{-10}$ and then to $8.10\times10^{-18}$, accepting a full step both times. The last residual is already below the $1.15\times10^{-15}$ stopping threshold. The exact velocity and pressure lie in the P2/P1 trial spaces, and the computed fields recover them to near floating-point precision. The Jacobian action agrees with centered differences to $8.6\times10^{-15}$ relative at the largest tested perturbation. This supports the nonlinear assembly and Jacobian action for this case. The gap is the length of the usable residual history.

**Follow-up.** Keep the same manufactured case, mesh, and boundary data. Start several Newton runs from controlled perturbations of the converged discrete solution. Apply perturbations to the free state variables so the prescribed velocity remains exact, and preserve the pressure gauge. Record accepted step lengths and residuals. Estimate the local order from consecutive residuals on full steps that are above the numerical floor. Close the rung once the runs show an order near two before roundoff dominates, using an acceptance band agreed before running the test.

## Rung 6. Cross-check the tangent with an adjoint

After converting the reported derivatives to the hand reference's dimensionless convention, the flow sensitivities agree within $2.0\times10^{-7}$ absolute. The power sensitivities agree as well. For depth, a tenfold reduction in perturbation reduces the corrected Taylor remainder by about 100, as expected for a second-order remainder. This supports the tangent implementation for flow and power on the active-drag strip. The existing runs compare tangents with centered differences and the continuum reference. An adjoint calculation is the missing independent check.

**Follow-up.** On the same converged strip and mesh, solve the transpose of the reduced state Jacobian once for each response, flow and power. Compute each total derivative using the residual parameter derivative and the direct dependence of the response on that parameter. This includes the explicit depth factor in flow and the driving-force and length factors in power. Compare adjoint, tangent, and centered differences in the same units, and repeat the Taylor checks. Use an absolute tolerance for derivatives whose reference is zero or near zero, such as flow with respect to strip length at fixed driving force. Record the tolerances before running these checks. Wall-shear sensitivities and the optional fixed-flow variant are separate follow-up coverage.

## What this shows about the solver

The tests support the verification solver's strip flow and power predictions, the expected convergence of the P2/P1 discretization, rotation consistency, and the nonlinear and tangent calculations exercised here. Wall shear is the least accurate strip quantity measured, and its error is consistent with the observed force imbalance because both measurements use the same reconstructed wall stress. The evidence covers the periodic strip and manufactured square. The production manifold and its outlet boundary conditions need their own verification.


## Follow-up implementation and results for Rungs 5 and 6

**Rung 5 PASS. Rung 6 PASS.** The verdicts above document the original review. This section records the first follow-up run.

The follow-up uses the original manufactured square and active-drag strip. The Newton routine in `fem.py` now accepts a starting reduced state. The runner uses that option for the one-step study, assembles the flow and power response derivatives, and calls the existing adjoint solver. Comments in both files explain the assembly, boundary reduction, measurements, and checks.

### Rung 5. Local Newton order

The study starts from the converged discrete state plus a reproducible perturbation in the free velocity unknowns. The starting perturbation preserves the prescribed boundary velocities and pressure gauge. The error is the RMS of the recovered velocity and pressure coefficient errors, scaled by the manufactured velocity and pressure scales. The study measures the remaining error after one full Newton step as the starting perturbation is halved.

The settings were fixed before running. Each direction uses seed 317 or 911 and amplitudes from 0.5 to 0.015625 in successive halvings. To pass, each seed must show an observed order between 1.8 and 2.2 on at least two consecutive intervals with full steps and final errors above the dimensionless floor of $10^{-10}$.

| Seed | Full steps | Usable observed orders | Verdict |
| --- | --- | --- | --- |
| 317 | 6/6 | 1.98527, 1.99553, 1.99855, 1.99947, 1.99979 | PASS |
| 911 | 6/6 | 1.99573, 1.99743, 1.99873, 1.99938, 1.99970 | PASS |

The measured orders establish local quadratic behavior for this manufactured case.

### Rung 6. Adjoint comparison

The runner solved the transpose of the reduced Jacobian once for flow and once for power. It included both the residual sensitivity and the direct dependence of each response on the parameter. Comparisons cover all nine strip inputs and use the absolute difference in derivatives normalized by the baseline parameter and response scales. This scale also applies to derivatives expected to be zero.

| Check | Measured result | Criterion set before the run |
| --- | --- | --- |
| Maximum tangent–adjoint difference | 2.334e-15 | At most 1e-8 |
| Maximum adjoint–centered difference at relative step 1e-4 | 5.595e-09 | At most 1e-5 |
| Depth Taylor orders, both responses and perturbation signs | 1.99910 to 2.00089 | All eight orders from 1.8 to 2.2, using steps 1e-2, 1e-3, and 1e-4 |

These checks cover flow and power at the fixed-driving-force baseline. Wall-shear sensitivities and the optional fixed-flow case remain outside this follow-up.

### Regression comparison

Fresh runs of the original and updated runner gave identical values for all previously reported results, including the four earlier rungs, the standard nonlinear solve, and the tangent and finite-difference checks.

The [raw follow-up results](verification_followup_results.json) contain the starting and final errors, Newton step histories, and per-parameter derivative comparisons. The verdicts apply to the verification cases.
