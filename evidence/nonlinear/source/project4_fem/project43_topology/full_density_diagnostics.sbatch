#!/bin/bash -l
#SBATCH --job-name=p43-full-density
#SBATCH --partition=htc
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --array=0-6%4
#SBATCH --output=p43-full-density-%A_%a.log

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
root="${P43_STUDY_ROOT:-/scratch/$USER/p43-topology-study-20261125}"
deps="${P43_DEPS:-$root/gate2-source-2395294/dependencies}"
python_path="${P43_PYTHON:-$HOME/.conda/envs/brinkman-rans/bin/python}"
source_dir="$root/full-density-source-$SLURM_JOB_ID"
results="$root/full-density-diagnostics-$SLURM_ARRAY_JOB_ID"
mkdir "$source_dir"
mkdir -p "$results"
tar -xzf "$SLURM_SUBMIT_DIR/full-density-diagnostics.tar.gz" -C "$source_dir"
cd "$source_dir"
export PYTHONPATH="$deps${PYTHONPATH:+:$PYTHONPATH}"
run() { /usr/bin/time -v "$python_path" -B -u -m "project4_fem.project43_topology.$@"; }

case "$SLURM_ARRAY_TASK_ID" in
  0) run connectivity "$results/connectivity.json" --archive "$root/recovery-results-2396219/mma/mma_summary.json" ;;
  1) run corrected_brinkman_channel "$results/corrected_channel.json" ;;
  2) run expanded_envelope "$results/envelope_baselines.json" ;;
  3) run full_density_continuation "$results/original_projection" --route projection --mma-evaluations 150 ;;
  4) run full_density_continuation "$results/original_brinkman" --route brinkman --mma-evaluations 150 ;;
  5) run full_density_continuation "$results/expanded_projection" --route projection --envelope expanded --mma-evaluations 150 ;;
  6) run full_density_continuation "$results/original_coupled" --route coupled --mma-evaluations 150 ;;
esac
