# Coolant manifold research campaign handoff

**Finding.** Moving half of the existing shoulder-pocket compensation downward reduced Re 50 side flow from **21.2062% to 21.1541%** on the fine fitted mesh. That is a repeatable local gain of **0.0520 percentage points**, leaving **0.6541 points** to reach the upper flow limit. Fine-grid pressure loss was **95.2946**, just **0.0252** below the provisional cap. The same design missed the cap on the medium grid. No tested point met the flow band. The results support a local tradeoff and leave feasibility across the full geometry family open.

## Experiment and reproducibility

The question was whether the best existing metering design could direct less coolant to the side outlet within its pressure margin. Four diagnostic coordinates shortened the plateau, shifted a fraction of the area compensation between the two existing pockets, and independently blended cubic and quintic entrance and exit ramps. The nominal width stayed at 1.00 mm and the design stayed within the existing supports with the same ports, boundary conditions and area rule. The pressure screen remained **95.31983921402909** and the side-flow band **19.5%–20.5%**. The independent compliance review approved one bounded stencil and a conditional held-out evaluation before submission. [Campaign scope](RESEARCH_CAMPAIGN_SCOPE.md) · [Frozen protocol](research_geometry/PROTOCOL.md) · [Independent compliance record](research_compliance/REVIEW.md)

Sol job **2085462** used the existing Stokes-initialized Re 50 Taylor–Hood fitted-wall solver. All **11 new L2** and **3 new L3** states converged. The maximum reported free-DOF residual was **8.46×10⁻¹⁴** and relative mass imbalance **1.14×10⁻¹⁴**. All converged shapes passed the width, port and area checks. Two origin states were reused from job 2084555 after verifying source hashes, the reference parameterization, and identical mesh connectivity and coordinates within the declared tolerance. The solver did not change. [Raw case rows, prediction errors and mesh fingerprints](research_geometry/sol-2085462/results/local-stencil.json)

The job completed in **100 s** with four workers, using **0.1111 allocated CPU-hours** and **0.06676 measured CPU-hours**. Slurm recorded **4.03 GiB** peak batch RSS and zero eligibility-to-start delay at one-second resolution. Individual L2 solves took **10.48–11.16 s**, about a quarter of the **39.51–40.54 s** L3 solves. The frozen protocol understated the L2 time relative to L3. The actual job was within its planned budget. The immutable source archive includes the exact solver dependencies, reference data, reviewed protocol, and executable configuration. The Sol environment was Python **3.12.16**, NumPy **2.5.4** and SciPy **1.18.1**. The base Git revision was **917c4600f7d95b28a19257110eea59e16d55e1ca**. [Source and exact inputs](research_geometry/sol-2085462/source.tar.gz) · [SHA-256 manifest](research_geometry/sol-2085462/input-sha256.txt) · [Independent HPC accounting](research_hpc/LOCAL_METERING_ACCOUNTING.md)

## Measured geometry and pressure tradeoff

Positive headroom means the measured pressure is below the cap. Paired side-flow changes use the origin on the same mesh. Pressure in the table is nondimensional.

| Fine-grid case | Side flow | Side-flow change (pp) | Pressure loss | Pressure headroom |
| --- | ---: | ---: | ---: | ---: |
| Upper-cubic origin | 21.2062% | — | 95.2336 | +0.08622 |
| Mixed held-out direction | 21.2851% | +0.07891 | 95.1247 | +0.19516 |
| **50% transfer to the lower pocket** | **21.1541%** | **−0.05201** | **95.2946** | **+0.02522** |
| 80% of the affine prediction | 21.1640% | −0.04215 | 95.2831 | +0.03675 |

The largest sampled single-axis reduction came from moving half the compensation from the upper to the lower shoulder pocket. The fine polygonal area gain was **+0.000978 mm²**, with a reconstructed gain of **+0.001000 mm²** relative to the matched 26-control baseline. The geometry kept its original 1.50 mm plateau and both cubic ramps. The two independent ramp changes each reduced side flow by about **0.017 points** at full L2 amplitude and each raised pressure by about **0.021–0.022**. Their responses have the same sign, so the tested direction gives no evidence of cancellation. Shortening the plateau lowered pressure but increased side flow. Every measured point with lower side flow than its same-grid origin paid a pressure increase. [Geometry and sensitivity details](research_geometry/RESULTS.md)

All six signed response comparisons for the three refined directions met the declared two-grid resolution rule. The largest changes in a paired response from L2 to L3 were **0.000328 points** for side flow and **0.000441** for pressure. Those small paired variations do not guarantee the absolute pressure classification. For the 50% pocket transfer, absolute pressure changed by **0.02655** between grids, larger than its fine headroom. It was **0.00133 over** the cap at L2.

The independent mixed point was **0.00030 points higher** in side flow and **0.00134 lower** in pressure than its L2 affine prediction. The locally predicted cap-bound point lowered side flow to **21.1619%** at L2 but exceeded pressure by **0.00412**, so the preregistered 80% step was evaluated and refined. The 80% step measured **21.1640% / 95.2831** on L3. The finite secants and this predictor resolve a local trend. They do not establish fitted-adjoint correctness or optimization stationarity. [Figure of measured flow and pressure](research_geometry/sol-2085462/results/local-authority.png) · [Independent numerical audit](research_verification/VERIFY.md)

## Independent pressure and RANS evidence

Archived upper-cubic and lower-quintic L2 cases were reanalyzed without repeating the successful OpenFOAM comparison. OpenFOAM pressure exceeded fitted FEM by **0.21531** and **0.40679**, respectively. For the same geometry change, FEM added **0.23127** pressure units and OpenFOAM added **0.42275**. The **0.19147** difference in these increments is evidence against a fixed solver offset. The outlet formulations and the discretizations both differ. The archive cannot isolate their contributions or give an error bound for a new shape. Both archived codes gave side flow above the target for both shapes. The smallest new FEM target gap, **0.6541 points**, exceeds the largest archived metering cross-code side-flow discrepancy, **0.1272 points**. An independent pressure reclassification alone would therefore leave every measured candidate outside the flow band. The protocol did not trigger a fresh OpenFOAM calculation. [Targeted OpenFOAM audit](research_openfoam/PRESSURE_AUDIT.md)

The published Project 4 manifold example allows Stokes or Navier–Stokes with Brinkman and names inertial nonlinearity in the Re 10–200 range. The Re 50 work supports that laminar route. The prior manufactured RANS residual and Jacobian verification passed, while the original adjoint result passed **10 of 12** checks. Smaller diagnostic steps improved agreement where limiter branches affected the finite differences. The original gate remains failed. The separately driven fitted RANS channel converged with revised near-wall treatment but had a **2.32%** medium-to-fine bulk change against the **1%** threshold. These results do not verify a physical turbulent-manifold adjoint. The RANS work remains conditional on an actual turbulent claim or a separate course requirement. [Official Project 4 specification](https://designinformaticslab.github.io/DesignOptimization2025/project4.html) · [Prior RANS evidence](next_rans_adjoint/README.md) · [Independent verification audit](research_verification/VERIFY.md)

## Research decisions

| Priority | Decision |
| --- | --- |
| **Required before a stronger design claim** | Establish the hardware head allowance and pressure measurement definition, available shoulder stock and cutter clearance, dimensional width tolerance, and the production area rule. The current pilot values can screen numerical designs. |
| **Conditional geometry investigation** | If the flow band remains mandatory, authorize a physically justified new source of geometric control or a wider tested region. The local stencil shows measurable but small gains as pressure headroom shrinks. A further sweep of the same local box lacks a distinct scientific question. |
| **Conditional OpenFOAM pressure investigation** | When a fresh fine-grid candidate reaches the flow band or a pressure classification becomes decision-critical, run a targeted cross-check with the accepted mesh family. Distinguish outlet-boundary effects from mesh sensitivity by controlling one change at a time. |
| **Required for eventual optimization evidence** | Verify the gradient for the actual selected geometry representation and all active responses, obtain a suitably specified feasible problem, and establish constrained stationarity and the required comparison with and without the nonlinear component. This local stencil adds design evidence without fulfilling that course deliverable. |
| **Conditional RANS campaign** | If required, establish a separate physical turbulent regime, resolve wall treatment and channel refinement, and pass the coupled derivative gates on representative physical states before using a RANS adjoint for design. |

Geometry and numerical verification were dispatched as GPT‑6.2 Astra with ultra requested. OpenFOAM and independent compliance were dispatched as GPT‑6.2 Astra with high requested. HPC inherited its launch settings, and this lead session exposed no verified model or effort readback. The active-agent limit was handled by starting the OpenFOAM and HPC reviews as slots opened. No reviewers submitted redundant CFD.

To repeat the bounded study on Sol, extract the archived sources into a fresh directory and submit:

```bash
sbatch project4_fem/research_hpc/run_local.sbatch project4_fem.research_geometry.local_stencil --output results --workers 4
```
