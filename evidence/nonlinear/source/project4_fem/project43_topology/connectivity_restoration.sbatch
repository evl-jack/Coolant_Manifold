#!/bin/bash -l
#SBATCH --job-name=p43-connectivity-restoration
#SBATCH --partition=htc
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --array=0-1
#SBATCH --output=p43-connectivity-restoration-%A_%a.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
root="/scratch/$USER/p43-topology-study-20261125"
source_dir="$root/connectivity-restoration-source-$SLURM_JOB_ID"
results="$root/connectivity-restoration-$SLURM_ARRAY_JOB_ID"
mkdir "$source_dir"
mkdir -p "$results"
tar -xzf "$SLURM_SUBMIT_DIR/connectivity-restoration.tar.gz" -C "$source_dir"
cd "$source_dir"
export PYTHONPATH="$root/gate2-source-2395294/dependencies"
envelopes=(original expanded)
envelope="${envelopes[$SLURM_ARRAY_TASK_ID]}"
"$HOME/.conda/envs/brinkman-rans/bin/python" -B -u -m project4_fem.project43_topology.full_density_continuation \
  "$results/$envelope" --route projection --envelope "$envelope" --mma-evaluations 150 --connectivity-diagnostic