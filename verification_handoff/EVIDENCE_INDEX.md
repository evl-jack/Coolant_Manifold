# Evidence index
Pinned commit bba3589e359daa6489b13446155b677bffd31b91. Core FEM shared Git blob 2826e76dc590d8f2d6d23c09f52a6f684cdffc4f.
| Evidence | Path | Job |
|---|---|---|
|Historical ladder|historical/baseline_verification.json; evidence/verification/source/project4_fem/verify.py|unknown|
|Pilot config/summary|evidence/verification/input/configuration.json; evidence/verification/input/summary.json|local|
|Input fingerprints|evidence/verification/archive_fingerprints.json|28 common inputs|
|Mapping|evidence/verification/source/project4_fem/project43_topology/fixed_pipe_pilot_model.py; evidence/verification/source/project4_fem/project43_topology/wall_displacement_density_v1.py; evidence/verification/source/project4_fem/project43_topology/density.py; evidence/verification/source/project4_fem/project43_topology/evaluator.py|local|
|Raw fields|evidence/verification/input/stage_03_beta_8/original/fields.npz; evidence/verification/input/stage_03_beta_8/selected/fields.npz|local|
|Raw contour|evidence/verification/input/stage_03_beta_8/selected/geometry.json|local|
|Gradient|evidence/verification/adjoint-evidence/adjoint_merged.json; evidence/verification/adjoint-results/|2408528,2408552|
|Supplement|evidence/verification/adjoint-coordinate-supplement.json|2408595,2408613|
|Timing|evidence/verification/timing-2408534/|2408534|
|KKT|evidence/verification/kkt-results/fixed_pipe_kkt.json|2408529|
|Mesh|evidence/verification/mesh-evidence/mesh_report.json; evidence/verification/mesh-results/|2408522 parent; individual IDs in JSON;2408535 report|
|Nonlinear|evidence/nonlinear/nonlinear_component.md; evidence/nonlinear/nonlinear_comparison.json; evidence/nonlinear/cross_evaluation/cross_evaluation.json; evidence/nonlinear/sol_refinement/stokes_slurm_accounting.txt|accounting lists jobs|
|Historical shapes|evidence/verification/source/project4_fem/MATHEMATICAL_OPTIMIZATION_HANDOFF.md; evidence/verification/source/project4_fem/TRANSITION_PILOT_HANDOFF.md|2085839,2085871|
|Historical topology|evidence/verification/source/project4_fem/project43_topology/gate2_evidence.md; evidence/verification/source/project4_fem/project43_topology/gate3_evidence.md|2395055–2396361|
|Comparator|evidence/verification/source/project4_fem/project43_topology/conventional_baseline.py; evidence/verification/source/project4_fem/project43_topology/wall_displacement_conventional_v1.py|not matched|
|Historical draft|historical/Project4_Report_Draft.md|superseded|
|Project uploads|Project43_searchable.txt, Project43_Sol_searchable.txt, Project43_SST_Newton_searchable.txt, report.json|headings and jobs in audit|
See [figure list](WEBSITE_VERIFICATION_HANDOFF.md), [raw tables](VERIFICATION_TABLES.json). Compressed uploads/uncommitted Sol scratch not exhaustively inventoried.
