#!/bin/bash -l
#SBATCH --job-name=p43-topology-restoration
#SBATCH --partition=htc
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=p43-topology-restoration-%j.log

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
study_root="${P43_STUDY_ROOT:-/scratch/$USER/p43-topology-study-20261125}"
dependency_dir="${P43_DEPS:-$study_root/gate2-source-2395294/dependencies}"
export PYTHONPATH="$dependency_dir${PYTHONPATH:+:$PYTHONPATH}"
python_path="${P43_PYTHON:-$HOME/.conda/envs/brinkman-rans/bin/python}"
job_source="$study_root/restoration-source-$SLURM_JOB_ID"
output="$study_root/restoration-results-$SLURM_JOB_ID"
mkdir "$job_source" "$output"
tar -xzf "$SLURM_SUBMIT_DIR/p43-topology-restoration.tar.gz" -C "$job_source"
cd "$job_source"
for mode in analytical full reduced; do
  /usr/bin/time -v "$python_path" -B -u -m project4_fem.project43_topology.gate2_restoration \
    "$output/$mode" "$mode"
done
