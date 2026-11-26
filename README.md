# Project 4 — GitHub Pages ready site

This is the finished static report and viewer implementation, with Markdown source and frozen numerical evidence. No build step is needed for publication.

## Import and publish
1. Extract this ZIP. Copy its **contents**, including `.nojekyll`, to the root of the public project repository. `index.md` is the main report source. Inspect the included evidence for any information the team does not want made public (historical host names and job paths are present).
2. Commit `index.html`, `index.md`, `assets/`, `data/`, evidence and source files. In GitHub Settings → Pages, select Deploy from a branch, main (or submission branch), `/ (root)`.
3. Open the deployed URL in a private browser. Check both WebGL modes, pressure and velocity, element click, stage playback, phone width and MathJax. Browser acceptance was blocked locally by environment policy; do not mark it passed until actually tested.
4. Add team authors and repository/site URLs, complete the AI disclosure from personal knowledge, and submit both URLs to Canvas as required by the [assignment](https://designinformaticslab.github.io/DesignOptimization2025/project4.html).

## Preview
`python -m http.server 8000` from this directory, then visit `http://localhost:8000/`. No npm or Jekyll install. Native WebGL and pinned MathJax 3.2.2 CDN. The latter requires network access.

## Rebuild from frozen data
`python -m pip install numpy shapely scikit-image markdown2` then `python rebuild_site.py`. This regenerates the consolidated JSON, viewer geometry, and report HTML from the supplied evidence and the main report in `index.md`. It does **not** rerun the physical simulations or dynamically regenerate every prose number. Do not represent this as satisfying a full end-to-end CFD reproduction requirement. The numerical source, scripts and original input/checkpoint data are in `evidence/`. See the [verification README](evidence/verification/README.md), [nonlinear README](evidence/nonlinear/README.md), and archived Sol launchers for actual rerun instructions. No new CFD has been performed to make the website.

## What remains
Read `FINAL_ACCEPTANCE.md` and `MISSING_DATA.md`. The archived baseline benchmark outputs are now included in `historical/baseline_verification.json`; historical benchmark rerun against final source and independent student calculations have not been established. The final model has one independent PDE inequality; no second constraint is claimed. The raw NS design is the primary candidate; smoothed fails the prescribed mesh screen. The conceptual cylindrical fluid-volume surface is not a 3D solution. The numerical display is a 0.6 mm extrusion of the actual 0.25 mm fitted analysis mesh, with 0.0625 mm performance data.
