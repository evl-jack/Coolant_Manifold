#!/bin/bash -l
#SBATCH --job-name=p43-topology-endpoint-derivatives
#SBATCH --partition=htc
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=p43-topology-endpoint-derivatives-%j.log

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
study_root="${P43_STUDY_ROOT:-/scratch/$USER/p43-topology-study-20261125}"
dependency_dir="${P43_DEPS:-$study_root/gate2-source-2395294/dependencies}"
export PYTHONPATH="$dependency_dir${PYTHONPATH:+:$PYTHONPATH}"
python_path="${P43_PYTHON:-$HOME/.conda/envs/brinkman-rans/bin/python}"
job_source="$study_root/endpoint-source-$SLURM_JOB_ID"
output="$study_root/endpoint-results-$SLURM_JOB_ID"
mkdir "$job_source"
tar -xzf "$SLURM_SUBMIT_DIR/p43-topology-endpoint.tar.gz" -C "$job_source"
cd "$job_source"
"$python_path" -B -u -m project4_fem.project43_topology.full_recovery \
  "$study_root/restoration-results-2395874/full" "$output" \
  --endpoint-root "$study_root/recovery-results-2396219"
