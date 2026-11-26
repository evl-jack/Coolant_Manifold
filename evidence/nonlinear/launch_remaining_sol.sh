#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
python="${P43_PYTHON:-/scratch/$USER/p43-accelerated/cpu-py312-jax053/bin/python}"
deps="${P43_DEPENDENCIES:-/scratch/$USER/p43-topology-study-20261125/gate2-source-2395294/dependencies}"
export PYTHONPATH="$root/source:$deps"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
stage="$root/archive/pilot/stage_03_beta_8"
contour="$root/stokes_run/stage_03_beta_8/selected/geometry.json"
outdir="$root/fitted_sol"
mkdir -p "$outdir"
job=$(sbatch --parsable --array=1-2 -o "$root/stokes-fitted-%A-%a.log" "$root/source/project4_fem/project43_topology/verification/stokes_fitted.sbatch" "$root" "$stage" "$contour" "$outdir")
cp -an "$root/fitted_stokes/fitted-stokes-0.25" "$outdir/"
report=$(sbatch --parsable -p htc -c 1 --mem=2G -t 00:10:00 -d "afterany:${job%%;*}" -o "$root/stokes-report-%j.log" --wrap "'$python' -B -m project4_fem.project43_topology.verification.stokes_fitted report '$stage' '$root/archive/fitted_mesh_results' '$outdir' '$root/fitted_sol_report'")
printf 'Fitted job %s and report job %s
' "$job" "$report"
