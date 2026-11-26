# Final grading checklist and remaining team actions

Reviewed official assignment: https://designinformaticslab.github.io/DesignOptimization2025/project4.html

| Category | Points | Included | Qualification / remaining action |
|---|---:|---|---|
| Formulation | 10 | explicit 128 controls, power, side flow, bounds, rationale | one independent PDE inequality; no invented second |
| Analysis and verification | 20 | weak form, mesh, exact fixture table, raw benchmark JSON, current verification | historical benchmark not rerun on final executable; confirm independently derived expected values |
| Sensitivities | 20 | both response derivatives, FD/Taylor checks incl. original flags, timing, scaling | no new run; use supplied evidence |
| Optimization | 15 | MMA, scaling, continuation, KKT and multiplier interpretation | local immersed first-order audit only |
| Nonlinear | 15 | matched Stokes/NS optimization and nonlinear evaluation, mesh, fields, differences | no global optimum claim |
| Website/3D | 10 | both geometry modes, linked controls, sampled FE fields, click inspection, actual checkpoint animation, metrics | code and static asset checks only; deployed browser, rendering, load time and phone acceptance pending |
| Reproduction | 5 | source, archived inputs/results, launch scripts, one-command archive-site rebuild | full end-to-end numerical regeneration not tested; some narrative values frozen |
| Clarity | 5 | full report, real figures, slide | final presentation and public deployment team action |
| N2 bonus | +5 possible | nonlinear benchmark and final nonlinear adjoints | instructor determines applicability of historical benchmark |
| Topology bonus | +5 possible | density filter, projection, topology change | equivalent fixed-topology comparator missing |

Eight viewer features are implemented in source: linked cameras, true extruded geometry, physical color scales, P2/P1 sampled fields, actual checkpoint playback, click inspection, separate metrics/status and JSON data. Native WebGL rendering and performance must be verified in a normal browser. The cylindrical surface is a sampled equivalent-width union, not an exact asymmetric contour or a 3D flow model. The playback is an actual checkpoint image sequence, not a time-dependent CFD animation.

No team names, human derivation claims, verified graded solver bug, second constraint, matched-area comparator, deployed URL, or actual awarded grade are invented.

## Local checks performed
- Python scripts executed to produce the report, viewer meshes and numerical JSON; syntax compilation passed.
- Static HTML local links verified (no missing linked files); 8 report tables and 19 image elements present.
- JavaScript syntax check passed. Independently checked P2 velocity-node correspondence with triangle vertices/midpoints for all four display meshes; no ambiguous element mapping.
- Browser navigation to the local preview was blocked by administrator policy. Live rendering, interaction, MathJax and phone layout are **not yet accepted**.
