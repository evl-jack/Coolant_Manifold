#!/bin/bash -l
#SBATCH --job-name=p43-fixed-adjoint
#SBATCH --partition=htc
#SBATCH --array=0-4
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

set -euo pipefail
source_dir="${1:?Pass the directory containing project4_fem}"
snapshots="${2:?Pass the pilot stage snapshot directory}"
output="${3:?Pass the audit output directory}"
python="${PILOT_PYTHON:-/scratch/$USER/p43-accelerated/cpu-py312-jax053/bin/python}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
export PYTHONPATH="$source_dir${PYTHONPATH:+:$PYTHONPATH}"
cases=(original gray beta_0 beta_2 beta_8)
case_name="${cases[$SLURM_ARRAY_TASK_ID]}"
mkdir -p "$output"
cd "$output"
exec "$python" -B -u -m project4_fem.project43_topology.verification.fixed_pipe_adjoint_audit \
    "$snapshots" "$output/adjoint-$case_name.json" --cases "$case_name"
