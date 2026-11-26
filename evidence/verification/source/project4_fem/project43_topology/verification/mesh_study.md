# Fixed-pipe mesh verification

Use the final archived `stage_03_beta_8` directory. Keep its `original` and `selected` snapshots intact. Supply the frozen `smoothed_geometry.json` exported by the smoothing study. Each run writes a SHA-256 digest of its inputs, and the report checks that these digests agree across refinement levels. The study uses the existing Taylor–Hood FEM, nonlinear residual, inlet profile, zero-total-traction outlets, port-power functional, and Triangle-based fitted mesher.

## Run on Sol

Stage the project source and inputs in this layout. Keep the parameterized mesh builder from `fixed_pipe_pilot_model.py` with the verification modules.

```text
$run_root/source/project4_fem/
$run_root/input/stage_03_beta_8/original/
$run_root/input/stage_03_beta_8/selected/
$run_root/input/smoothed_geometry.json
```

Submit the independent cases as one array. The output argument is an absolute path to a new audit directory. The defaults in the batch script match the existing Sol project environment. Set `P43_PYTHON` and `P43_DEPENDENCIES` if installations are staged elsewhere.

```bash
cd "$run_root"
mesh_job=$(sbatch --parsable --output="$run_root/mesh-%A-%a.log" \
  source/project4_fem/project43_topology/verification/mesh_cases.sbatch \
  "$run_root" "$run_root/mesh-results")
sacct -j "$mesh_job" --format=JobID,State,ExitCode,NodeList,Elapsed,TotalCPU,MaxRSS
```

The array replays the original and selected immersed fields and the original, raw, and smoothed physical contours at analysis targets of 0.25, 0.125, and 0.0625 mm. The immersed design grid stays at 0.25 mm. The fitted geometry remains the same polygon at every level. Junction and collar maximum triangle edges quantify local resolution in addition to the nominal mesh target. An optional fourth level at 0.03125 mm is available with `--array=15-19%5`, contingent on the first report. Keep the initial report and task logs if adding this level.

After the array finishes, regenerate all mesh tables and plots with one command. Supply `--additional-level` after running the fourth level.

```bash
dependencies="${P43_DEPENDENCIES:-/scratch/$USER/p43-topology-study-20261125/gate2-source-2395294/dependencies}"
python="${P43_PYTHON:-/scratch/$USER/p43-accelerated/cpu-py312-jax053/bin/python}"
export PYTHONPATH="$run_root/source:$dependencies" MPLBACKEND=Agg
"$python" -B -m project4_fem.project43_topology.verification.mesh_report \
  "$run_root/mesh-results" "$run_root/mesh-evidence"
```

The case table and JSON include failed and missing cases. Failed solves retain their diagnostic and exit nonzero. Figure field maps use the finest successful result for each geometry. The color scales are common within the immersed and fitted families, and the velocity plot subdivides each P2 element at its midpoint nodes for rendering.

## Energy and leakage conventions

Port power is the signed pressure and kinetic energy entering across all three ports. The raw energy residual sums that power and the directly integrated viscous stress work at the ports, then subtracts viscous and Brinkman dissipation. The divergence-corrected result adds the pressure-plus-kinetic weighted divergence integral. The independent weak closure tests the momentum residual against the finite-element velocity, retaining reactions on prescribed boundary degrees of freedom. Pairwise comparisons report the sum of absolute raw residual magnitudes alongside the power saving.

The immersed transport diagnostic reports RMS flow and gross ingress and egress across interfaces of element-aligned projected-density cells below 0.5. It measures transport in the penetrable modeled solid. Its cell mask differs from the interpolated physical wall contour. The minimum nominal penetration length is 0.10010 mm from `sqrt(mu/alpha_max)`. The three immersed grids provide approximately 0.40, 0.80, and 1.60 nominal steps per minimum penetration length. Mesh convergence of the finite-resistance PDE and resolution of its solid penetration are separate questions.

## Decision rule

Each individual series must meet both successive-difference thresholds over its two latest refinement pairs. The thresholds are 0.5% for mechanical power, relative to the finer mesh, and 0.020 percentage points for side flow. The comparison tracks the signed paired power saving and side-flow change. Its uncertainty screen uses the maximum over the two latest absolute paired changes and the corresponding sums of individual absolute changes. This is descriptive numerical uncertainty. Physical savings require convergence screens to pass, a positive finest saving above the paired estimate, a combined raw energy residual smaller than the saving, and a side-flow increase large enough to cover the corresponding side-flow estimate. For a monotone sequence with decreasing increments, an apparent order and extrapolated value are reported as diagnostics for the nominal 2:1 refinement. Numerical differences across the available grids remain the basis of the decision.

The cases use the same inlet volume flux and produce different outlet shares and physical fluid areas. Estimating a curvature-only effect would require matched-area and matched-side-flow controls.
