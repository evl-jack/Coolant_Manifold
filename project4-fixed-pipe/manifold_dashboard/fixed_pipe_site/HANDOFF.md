# Project 4 website handoff

## Build and preview

From the package or project root with Python 3.12 or later, install the build dependencies and generate the report.

    python -m pip install -r manifold_dashboard/fixed_pipe_site/requirements-build.txt
    python manifold_dashboard/fixed_pipe_site/build.py
    python -m http.server 8000 --directory manifold_dashboard/fixed_pipe_site/dist

Open [http://localhost:8000/](http://localhost:8000/). The build produces HTML, MathML, field and geometry assets, playback data and the numerical manifest. WebGL and browser-native MathML require no runtime dependencies from a CDN.

The packaged **dist/offline-preview.html** also opens directly from disk. Keep it alongside **dist/assets** for images and links. Run the build with **--offline-preview** to regenerate it.

## Sources used

| Source | Use |
|---|---|
| Official ASU Project 4 assignment | Report, viewer, rubric and submission requirements |
| Pasted text(8).txt | Fixed-pipe model and deliverable scope |
| output/fixed_pipe_pilot_local/configuration.json and summary.json | Geometry bounds, Reynolds number, design count, MMA stages |
| output/fixed_pipe_pilot_local/mesh.npz and stage geometry and field archives | Numerical extrusion, branch widths, original and selected Taylor–Hood states, checkpoint playback |
| output/fixed_pipe_pilot_local/progress.jsonl and metrics.json | MMA history, objective, side-flow limit and performance |
| output/fixed_pipe_fitted_smoothing/results.json and three_panel.png | Preliminary fitted-wall comparisons |
| output/fixed_pipe_kkt_endpoint/fixed_pipe_kkt.json | Local β = 8 endpoint replay |
| output/fixed_pipe_independent_verification/adjoint_taylor_remainders.png | Local gradient smoke fixture |
| Project4_Report_Draft.md | Common FEM verification and reported AI-use evidence |
| project4_fem/fem.py, junction.py, project43_power.py and project43_topology modules | Interpolation, state indexing, flux and power formulas, density pipeline |

Downloads for the linked Fixed_Pipe_Pilot.zip and Project4_Report_Package.zip could not complete. The build uses the pilot files already in the workspace and the report draft read directly. The workspace comparison PNG matches the linked filename and size. The earlier dashboard uses a different parameterization.

## Insert released numerical evidence

Place a reviewed JSON object at **output/fixed_pipe_site_release.json**, or in the directory passed to **--evidence**. Each key below maps to a record with **source**, **experiment**, **status**, **settings**, **columns** and **rows**. Put units in the column headings and use null for unavailable numbers. Table rows must match the column count. Optional **figure** is a raster filename beside the release JSON. **figureAlt** describes that figure.

| Key | Expected numerical source |
|---|---|
| finalTaylor | Both response gradients, finite differences and Taylor orders for feasible perturbations |
| cost | Matched hardware, assembly, factorization, solve counts and wall time |
| scaling | State size, design dimension, wall time and memory by mesh |
| independentKkt | Endpoint ID, active set, scaled multipliers and residuals |
| immersedMesh | Mesh settings and paired power and branch-share changes |
| fittedMesh | Fixed contours, paired meshes, power and flow uncertainty |
| stokes | Independently optimized Stokes design and matched nonlinear cross-evaluation |

Run the build again and review the tables, figure links and statuses. Retain the original numerical evidence and its experiment identifiers.

## Numerical reproduction

The site build reproduces the presentation from archived inputs. The pilot optimization can be rerun into a new output directory using the supplied solver source and compatible Python packages.

    python -m pip install -r manifold_dashboard/fixed_pipe_site/requirements-numerical.txt
    python -m project4_fem.project43_topology.fixed_pipe_pilot output/fixed_pipe_reproduction --local

The archived smoothing driver reads and writes fixed relative paths. Use a separate working copy if repeating the fitted calculation so that the supplied result remains intact.

Complete numerical reproduction additionally needs the paired fitted-wall, mesh and independent verification release inputs.

## GitHub Pages and Canvas

1. Use a clean public submission repository and review the source and evidence to be published, along with team attribution.
2. Copy **pages-workflow.yml** into **.github/workflows/pages.yml**. In the repository, select **Settings → Pages → GitHub Actions**.
3. Commit and push the reviewed package, then run **Project 4 Pages** for the grading revision.
4. Open the published URL in a private window. Check it at desktop and phone widths, including viewer controls, numerical assets, equations and links.
5. Submit the public repository and Pages links on Canvas. Retain the grading commit hash.

## Assignment checklist

| Requirement | Status |
|---|---|
| Headline results and comparison at top | Implemented from preliminary fitted-wall data |
| Motivation, formula, bounds and source of limits | Implemented for the pilot, which has one PDE-based inequality |
| Governing PDE, FEM, mesh, loads and assumptions | Implemented |
| At least six verification rungs | Historical common FEM measurements included |
| Adjoint, finite differences and Taylor curves | Derivations and local fixture included. Independent final evidence pending |
| Measured adjoint cost and scaling | Awaiting released data |
| MMA history, active constraints and multipliers | Archived history and local endpoint replay included. Independent release pending |
| Stokes versus nonlinear optimum | Awaiting fixed-pipe Stokes release |
| Objective and all constraints in results table | Immersed and fitted results included with distinct labels |
| Bonus evidence | N2 and topology scope documented |
| Limitations and AI-use disclosure | Draft included. Team attribution and main solver regression example pending |
| Linked cameras and true extrusion | Implemented, including 0.8 mm display depth |
| Cylindrical engineering interpretation | Implemented from measured transition widths |
| Fields, color scales and element interpolation | Implemented for β = 8 immersed velocity and pressure |
| Animation and inspection | Saved checkpoints and element hover/click implemented |
| Objective and constraint pass/fail | Implemented for both views |
| JSON, mobile responsiveness and loading | Responsive views and lazy field loading |
| Markdown, equations and assignment link | Included, math compiled to native MathML |
| Source code and one command to regenerate results.json | Archived-input build and standalone pilot entry point. End-to-end CFD reproduction awaits evidence |
| Public repository, deployed site, private browser, Canvas and grading commit | Requires public release and course submission |
