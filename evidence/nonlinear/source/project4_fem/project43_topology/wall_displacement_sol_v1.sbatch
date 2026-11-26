#!/bin/bash -l
#SBATCH --job-name=p43-wall-displacement-v1
#SBATCH --partition=htc
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=wall-displacement-v1-%j.log

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
root="$P43_WALL_ROOT"
source="$root/source-v1"
mode="$P43_WALL_MODE"
output="$root/$mode-$SLURM_JOB_ID"
python="/scratch/$USER/p43-accelerated/cpu-py312-jax053/bin/python"
export PYTHONPATH="$source:/scratch/$USER/p43-topology-study-20261125/gate2-source-2395294/dependencies"
cd "$source"

case "$mode" in
  verify)
    /usr/bin/time -v "$python" -B -u -m project4_fem.project43_topology.wall_displacement_verify_v1 "$output"
    ;;
  conventional)
    /usr/bin/time -v "$python" -B -u -m project4_fem.project43_topology.wall_displacement_conventional_v1 "$output"
    ;;
  optimize)
    /usr/bin/time -v "$python" -B -u -m project4_fem.project43_topology.wall_displacement_optimize_v1 \
      "$output" --common "$P43_WALL_COMMON"
    ;;
  *)
    echo "Unknown wall-displacement mode: $mode" >&2
    exit 64
    ;;
esac
