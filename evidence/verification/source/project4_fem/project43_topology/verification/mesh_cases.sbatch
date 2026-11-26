#!/bin/bash -l
#SBATCH --job-name=p43-mesh-audit
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH --array=0-14%5

set -euo pipefail
root="${1:?Pass the staging directory}"
output="${2:?Pass a new audit output directory}"
python="${P43_PYTHON:-/scratch/$USER/p43-accelerated/cpu-py312-jax053/bin/python}"
dependencies="${P43_DEPENDENCIES:-/scratch/$USER/p43-topology-study-20261125/gate2-source-2395294/dependencies}"
export PYTHONPATH="$root/source:$dependencies"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg

cases=(immersed-original immersed-raw fitted-original fitted-raw fitted-smoothed)
spacing=(0.25 0.125 0.0625 0.03125)
case_name="${cases[SLURM_ARRAY_TASK_ID % 5]}"
h="${spacing[SLURM_ARRAY_TASK_ID / 5]}"
cd "$root"
"$python" -B -u -m project4_fem.project43_topology.verification.mesh_cases \
    "$root/input/stage_03_beta_8" "$root/input/smoothed_geometry.json" \
    "$output/$case_name-$h" "$case_name" "$h"
