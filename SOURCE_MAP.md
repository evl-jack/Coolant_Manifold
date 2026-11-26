# Authoritative source map
Paths are relative to this package. The original archive hashes and entry counts are in SOURCE_MANIFEST.json. All three archives are extracted without modifications. New website content must be derived from these, not old chat summaries.

| Section | Exact evidence path |
|---|---|
| Final result overview and mesh/adjoint/KKT qualification | `evidence/verification/README.md` |
| Nonlinear comparison narrative | `evidence/nonlinear/README.md`, `evidence/nonlinear/nonlinear_component.md` |
| Final design settings | `evidence/verification/input/configuration.json` |
| Pilot figures and histories | `evidence/pilot/local/stage_03_beta_8/images/`, `evidence/pilot/local/` |
| Optimizer history | `evidence/verification/input/progress.jsonl` |
| Original / primary physical geometry | `evidence/verification/input/stage_03_beta_8/{original,selected}/geometry.json` |
| Optional smoothed geometry | `evidence/verification/input/smoothed_geometry.json`, `assets/three_panel.png` |
| Adjoint derivation | `evidence/verification/source/project4_fem/project43_topology/verification/ADJOINT_METHOD.md` |
| Initial directional checks (retain failed Taylor flags) | `evidence/verification/adjoint-evidence/adjoint_merged.json`, `adjoint-evidence/adjoint_directional_derivatives.csv` and `adjoint-evidence/adjoint_taylor_remainders.png` under same root |
| Larger-step follow-up | `evidence/verification/adjoint-coordinate-supplement.json`, `evidence/verification/adjoint-coordinate-evidence/` |
| Complete gradients and timing | `evidence/verification/timing-2408534/{pilot-128,scale-32,scale-128,scale-512}/benchmark.json` |
| Timing figure | `evidence/verification/timing-2408534/figures/fixed_pipe_sensitivity_cost.png` |
| KKT | `evidence/verification/kkt-results/fixed_pipe_kkt.json` and PNGs in same directory |
| Mesh, energy, leakage | `evidence/verification/mesh-evidence/mesh_report.json`, `mesh_table.csv`, `mesh_convergence.png` under same directory |
| Fitted pressure/velocity figures | `evidence/verification/mesh-evidence/fitted_pressure_speed.png` |
| Stokes cross-evaluation and optimization | `evidence/nonlinear/nonlinear_comparison.json`, `evidence/nonlinear/cross_evaluation/` |
| Stokes fitted refinement | `evidence/nonlinear/fitted_report/fitted_comparison.json`, `fitted_comparison.csv` |
| Nonlinear figures | `evidence/nonlinear/figures/` (12 PNGs) |
| Stokes method | `evidence/nonlinear/source/project4_fem/project43_topology/verification/STOKES_METHOD.md` |
| Source implementations | `evidence/verification/source/` and `evidence/nonlinear/source/` (keep separate snapshots) |
| Reproduction | `evidence/pilot/README.md`, `evidence/verification/launch_sol_audit.sh`, `rebuild_figures.sh`, `evidence/nonlinear/launch_remaining_sol.sh` and README instructions |
| Historical mathematical draft and verification figures | `historical/Project4_Report_Draft.md`, `historical/Project4_Report_Package.zip` |
| Historical optimizer evidence | `historical/Project43_source_evidence_OPT678.zip` |

The historical report is outdated: its 416/176-variable geometry and 10% requirement are NOT the final study. Use it only for mathematical background, earlier tests after checking applicability, and design-history/failed-study discussion. The final archived configuration and solver code govern all numbers. The older report's benchmark values are not independently certified by this handoff; trace them to raw fixtures before marking them passed. Do not merge solver versions blindly. Keep provenance for each reported number and figure.
