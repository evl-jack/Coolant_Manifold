# Fixed-pipe adjoint and centered-difference costs

The timing study uses the archived β = 8 unchanged immersed reference for the power normalization and the side-flow requirement. It evaluates a smooth interior gray field sampled at the editable grid cell centers. The two responses are \(f=P/P_{\mathrm{ref}}\) and \(\eta=Q_{\mathrm{side}}/Q_{\mathrm{in}}\). The archived inequality derivative is \(-\nabla \eta/0.001\). Both responses are evaluated from each perturbed forward solution.

| Case | Density spacing (mm) | FEM spacing (mm) | Density controls | FEM triangles | FEM degrees of freedom | Centered FD solves including baseline |
|---|---:|---:|---:|---:|---:|---:|
| Pilot serial comparison | 0.25 | 0.25 | 128 | 2,848 | 13,359 | 257 |
| Fixed FEM scaling | 0.50 | 0.125 | 32 | 11,392 | 52,347 | 65 |
| Fixed FEM scaling | 0.25 | 0.125 | 128 | 11,392 | 52,347 | 257 |
| Fixed FEM scaling | 0.125 | 0.125 | 512 | 11,392 | 52,347 | 1,025 |

The three scaling cases use exactly the same FEM vertices and connectivity. All density faces coincide with FEM edges. The pilot comparison keeps the actual pilot FEM resolution. The archived optimization endpoint is never used as a perturbed design or changed by this study.

The common baseline starts cold and each perturbation starts from that converged gray-field state. The adjoint uses one forward state and one shared sparse factorization of the reduced tangent, followed by two separately timed transpose substitutions. The primary centered step is \(10^{-4}\). The pilot case also probes \(10^{-2}\) through \(10^{-6}\) on a reproducible coordinate and the largest-magnitude gradient coordinate for each response.

The launcher runs all cases on one Sol allocation with single-thread numerical libraries. The pilot comparison uses one process for each method. The scaling cases use independent FD workers. The report retains elapsed wall time and aggregate CPU and wall work, including worker launch overhead in parallel elapsed time. The one-worker pilot ratio is the measured serial speedup. Model construction and imports are excluded from the direct comparison. The JSON records model setup time, solver timings, factors, substitutions, filter and response work, gradient pullback work, Newton updates, residuals, and parent and worker peak RSS on Linux. Peak RSS is cumulative for each process. The FD and adjoint baseline timings share one measured nonlinear solve with separately measured response postprocessing.

## Sol execution

Stage the complete source under `$stage/source`. Set `$reference` to the archived `stage_03_beta_8/original/metrics.json` on Sol. Use the existing pilot Python environment or override it with `P43_PYTHON`.

```bash
sbatch --output="$stage/timing-%j.log" \
    "$stage/source/project4_fem/project43_topology/verification/fixed_pipe_timing.sbatch" \
    "$stage" "$reference"
```

The launcher writes `timing-<job-id>/` and produces the figure when every case succeeds. Each per-case `adjoint.json` and `full_fd_progress.jsonl` is available during execution. Each complete `benchmark.json` contains both full coordinate gradients and their errors. The supplemental `step_sweep.json` holds the perturbation-size series. The following command regenerates the figure.

```bash
cd "$stage/source"
python -m project4_fem.project43_topology.verification.fixed_pipe_timing \
    plot "$stage/timing-$job" "$stage/timing-$job/figures"
```

For a lightweight entry-point check, run `pilot-128` with `--coordinate-limit 2 --skip-sweep`. The results explicitly mark this as a partial gradient.
