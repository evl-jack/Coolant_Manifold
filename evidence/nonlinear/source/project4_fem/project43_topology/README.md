# Project 4.3 laminar topology study

This study compares a 20-control Wendland reduced-density model with one design variable per cell. Both models use the original sharp T-junction, the same fixed ports, a maximum fluid area of 83 mm², and a target of 49.9–50.1% side-outlet flow. Power is measured in watts per unit depth and normalized by a fixed common immersed baseline.

- `density.py` defines the two design maps and their shared Helmholtz, projection and Brinkman pipeline.
- `domain.py` and `evaluator.py` build the fixed union mesh and solve the planar FEM state and response adjoints.
- `restoration.py`, `optimize.py` and `kkt_audit.py` implement Phase I feasibility, Phase II power minimization, and the independent stationarity audit.
- `verify_density.py`, `verify_flow.py`, `kkt_fixtures.py` and `test_reconstruction.py` verify the numerical implementation. The Slurm scripts and their Python entry points separate the verification and physical screens.
- `reduced_authority.py` measures the attainable density ranges and tests the reduced design space. `full_recovery.py` replays independent endpoint audits and bounded full-space optimizer restarts.
- `reconstruction.py`, `leakage.py` and `local_refinement.py` implement the threshold geometry and the fixed binary fixture diagnostics.

Numerical experiments run on ASU Sol. Each submission preserves its source bundle and uses a fresh job-ID output directory. The two sequential approvals for production are the optimizer verification and the physical Brinkman mesh and resistance assessment.

| Record | Evidence |
|---|---|
| [Gate 1 review](gate1_review.md) | Approved formulation and verification requirements |
| [Gate 2 evidence](gate2_evidence.md) | Density, transfer, adjoint, restoration, reduced authority and optimizer results |
| [Restoration protocol](restoration_protocol.md) | Shared epigraph formulation and pilot budgets |
| [Gate 3 protocol](gate3_protocol.md) | Prospective fixtures, mesh changes, leakage thresholds and resolution requirements |
| [Gate 3 evidence](gate3_evidence.md) | Fitted references and the immersed resistance screens |

## Review after the bounded diagnostics

The [project review](https://chatgpt.com/g/g-p-5df5c39440c08191bf518d6456de7dbc-jacks-solver/c/6b05bfb1-9b00-83e8-8d32-af9075dbf801) accepted the tested density and adjoints, feasibility restoration, first-order MMA pilot, and fitted-wall references. SLSQP stationarity remains unresolved. The original reduced basis demonstrated insufficient flow authority at the pilot setting. The tested Brinkman configurations failed the physical acceptance gate, including the interface-refined binary fixtures. Production awaits a physically validated formulation.

The proposed follow-up preserves the original study and evaluates two distinct changes. The reduced model would retain 20 controls and test partition-of-unity normalization before any support redesign. The binary Brinkman benchmark would test a subcell displacement of the computational penalty mask against independent channel and fitted-manifold references. A differentiable extension and production-cost assessment would follow only if the binary correction meets the physical criteria.
