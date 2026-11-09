# Coolant Manifold

This repository contains the working files for a PDE-based design-optimization project on a coolant distribution manifold. The project is being developed in stages: first establish and verify the analysis, then add sensitivities and optimization, and later publish the report and design comparison as a GitHub Pages site.

## Project status

The current analysis plan is a fixed-geometry, steady, incompressible 2D depth-averaged flow model with one inlet and a configurable number of outlets. Outlet flow rates are measured against the equal-flow target, (Q_{\mathrm{in}}/n). Brinkman penalization and RANS are part of the broader project plan, but a production manifold solver, RANS implementation, and optimization loop are not yet complete.

The numerical verification folder records a six-rung ladder for the benchmark implementation. Its latest status and measurements are documented in [verification status](Manifold%20Model/numerical_verification/verification_status.md) and [follow-up results](Manifold%20Model/numerical_verification/verification_followup_results.json). These results verify the benchmark cases; they do not establish that a complete manifold-flow solver is finished.

## Repository contents

- [Problem statement](Outline/Problem.md) — project motivation and goals.
- [Analysis outline](Outline/Analysis.md) — governing equations, boundary conditions, and modeling assumptions.
- [Manifold Model README](Manifold%20Model/ReadMe.md) — implementation scope and next steps.
- `Manifold Model/core/` — geometry, material, boundary-condition, Brinkman-term, and implicit-sensitivity building blocks. Some dimension-specific element and RANS folders are placeholders.
- `Manifold Model/numerical_verification/` — verification code, independent hand and numerical writeups, status history, and recorded data.
- `Manifold Model/docs/` — hand-reference notes and research notes on 3D meshing and gradients.
- [Python dependencies](Manifold%20Model/requirements.txt) — NumPy and SciPy requirements.

The separate top-level `Solver/` folder from the working Project 4 directory is intentionally excluded because it is not the solver for this manifold model.

## Planned work

1. Complete and verify the fixed-geometry 2D manifold solver.
2. Establish and check sensitivities for every response and design parameter.
3. Add the optimization loop after the analysis and gradients are verified.
4. Develop the 3D model and assess fixed designs with an appropriately specified RANS model.
5. Build the GitHub Pages report and interactive 3D comparison.

The final course report will be curated in the MAE 494 repository later. The assignment's [report requirements](https://designinformaticslab.github.io/DesignOptimization2025/project4.html#what-the-report-must-contain) are the guide for organizing that report.
