#!/bin/bash -l
#SBATCH --job-name=p43-topology-local-binary
#SBATCH --partition=htc
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --array=0-1
#SBATCH --output=p43-topology-local-binary-%A_%a.log

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
study_root="${P43_STUDY_ROOT:-/scratch/$USER/p43-topology-study-20261125}"
dependency_dir="${P43_DEPS:-$study_root/gate2-source-2395294/dependencies}"
export PYTHONPATH="$dependency_dir${PYTHONPATH:+:$PYTHONPATH}"
python_path="${P43_PYTHON:-$HOME/.conda/envs/brinkman-rans/bin/python}"
job_source="$study_root/local-binary-source-$SLURM_JOB_ID"
output="$study_root/local-binary-results-$SLURM_ARRAY_JOB_ID"
mkdir "$job_source"
mkdir -p "$output"
tar -xzf "$SLURM_SUBMIT_DIR/p43-topology-recovery.tar.gz" -C "$job_source"
cd "$job_source"
fixtures=(sharp obstructed)
fixture="${fixtures[$SLURM_ARRAY_TASK_ID]}"
/usr/bin/time -v "$python_path" -B -u -m project4_fem.project43_topology.gate3_local \
  "$fixture" 1e7 "$study_root/gate3-fitted-2395603/fitted-0.0625.json" "$output/$fixture" \
  --base-spacing 0.125 --max-triangles 140000
