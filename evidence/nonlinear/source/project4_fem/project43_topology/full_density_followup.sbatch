#!/bin/bash -l
#SBATCH --job-name=p43-connected-corrected
#SBATCH --partition=htc
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=00:45:00
#SBATCH --array=0-5%4
#SBATCH --output=p43-connected-corrected-%A_%a.log

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
root="/scratch/$USER/p43-topology-study-20261125"
python_path="$HOME/.conda/envs/brinkman-rans/bin/python"
source_dir="$root/full-density-followup-source-$SLURM_JOB_ID"
results="$root/full-density-followup-$SLURM_ARRAY_JOB_ID"
mkdir "$source_dir"
mkdir -p "$results"
tar -xzf "$SLURM_SUBMIT_DIR/full-density-followup.tar.gz" -C "$source_dir"
cd "$source_dir"
export PYTHONPATH="$root/gate2-source-2395294/dependencies"
run() { /usr/bin/time -v "$python_path" -B -u -m "project4_fem.project43_topology.$@"; }
fitted="$root/gate3-fitted-2395603/fitted-0.0625.json"
case "$SLURM_ARRAY_TASK_ID" in
  0) run full_density_continuation "$results/original_connectivity" --route projection --mma-evaluations 150 --connectivity-diagnostic ;;
  1) run full_density_continuation "$results/expanded_connectivity" --route projection --envelope expanded --mma-evaluations 150 --connectivity-diagnostic ;;
  2) run corrected_brinkman_manifold sharp 1e6 "$fitted" "$results/corrected_sharp_1e6.json" ;;
  3) run corrected_brinkman_manifold obstructed 1e6 "$fitted" "$results/corrected_obstructed_1e6.json" ;;
  4) run corrected_brinkman_manifold sharp 1e7 "$fitted" "$results/corrected_sharp_1e7.json" ;;
  5) run corrected_brinkman_manifold obstructed 1e7 "$fitted" "$results/corrected_obstructed_1e7.json" ;;
esac