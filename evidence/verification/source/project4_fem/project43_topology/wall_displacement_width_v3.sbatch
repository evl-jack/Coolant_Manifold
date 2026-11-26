#!/bin/bash -l
#SBATCH --job-name=p43-wall-width-v3
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=wall-width-v3-%j.log

set -euo pipefail
mode="${1:?Pass verify or optimize}"
root="${2:?Pass the width-only experiment directory}"
checkpoint="${3:?Pass the verified Navier-Stokes beta=4 checkpoint}"
common="${4:?Pass the sharp fitted reference directory}"
shift 4

cd "$root/source"
export PYTHONPATH="$PWD:/scratch/$USER/p43-topology-study-20261125/gate2-source-2395294/dependencies"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
python="/scratch/$USER/p43-accelerated/cpu-py312-jax053/bin/python"
exec "$python" -B -u -m project4_fem.project43_topology.wall_displacement_width_v3 \
    "$mode" "$root/$mode-$SLURM_JOB_ID" --checkpoint "$checkpoint" --common "$common" "$@"
