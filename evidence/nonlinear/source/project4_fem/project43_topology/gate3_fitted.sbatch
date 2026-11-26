#!/bin/bash -l
#SBATCH --job-name=p43-topology-fitted
#SBATCH --partition=htc
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=p43-topology-fitted-%j.log

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLBACKEND=Agg
study_root="${P43_STUDY_ROOT:-/scratch/$USER/p43-topology-study-20261125}"
dependency_dir="${P43_DEPS:-$study_root/gate2-source-2395294/dependencies}"
export PYTHONPATH="$dependency_dir${PYTHONPATH:+:$PYTHONPATH}"
python_path="${P43_PYTHON:-$HOME/.conda/envs/brinkman-rans/bin/python}"
job_source="$study_root/gate3-source-$SLURM_JOB_ID"
output="$study_root/gate3-fitted-$SLURM_JOB_ID"
mkdir "$job_source" "$output"
tar -xzf "$study_root/gate3-fitted-input/p43-topology-gate3.tar.gz" -C "$job_source"
cd "$job_source"
for spacing in 0.25 0.125 0.0625; do
  /usr/bin/time -v "$python_path" -B -u -m project4_fem.project43_topology.gate3_screen \
    1e5 "$spacing" "$output/fitted-$spacing.json" --fitted-only
done
