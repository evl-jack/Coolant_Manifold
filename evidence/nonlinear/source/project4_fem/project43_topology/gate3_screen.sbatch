#!/bin/bash -l
#SBATCH --job-name=p43-topology-gate3
#SBATCH --partition=htc
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --array=0-8%3
#SBATCH --output=p43-topology-gate3-%A_%a.log

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLBACKEND=Agg
: "${P43_SOURCE:?Set P43_SOURCE to the unpacked reviewed study bundle}"
: "${P43_DEPS:?Set P43_DEPS to the isolated verified dependency directory}"
: "${P43_RESULTS:?Set P43_RESULTS to the screen output directory}"
: "${P43_FITTED_RESULTS:?Set P43_FITTED_RESULTS to the fitted comparison directory}"
export PYTHONPATH="$P43_DEPS${PYTHONPATH:+:$PYTHONPATH}"
python_path="${P43_PYTHON:-$HOME/.conda/envs/brinkman-rans/bin/python}"
penalties=(1e5 1e6 1e7)
spacings=(0.25 0.125 0.0625)
alpha="${penalties[$((SLURM_ARRAY_TASK_ID / 3))]}"
spacing="${spacings[$((SLURM_ARRAY_TASK_ID % 3))]}"
mkdir -p "$P43_RESULTS"
cd "$P43_SOURCE"
/usr/bin/time -v "$python_path" -B -u -m project4_fem.project43_topology.gate3_screen \
  "$alpha" "$spacing" "$P43_RESULTS/screen-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}.json" \
  --reference "$P43_FITTED_RESULTS/fitted-$spacing.json"
