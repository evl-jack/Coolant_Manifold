#!/bin/bash -l
#SBATCH --job-name=p43-topology-gate2
#SBATCH --partition=htc
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=p43-topology-gate2-%j.log

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLBACKEND=Agg
python_path="${P43_PYTHON:-$HOME/.conda/envs/brinkman-rans/bin/python}"
study_root="${P43_STUDY_ROOT:-/scratch/$USER/p43-topology-study-20261125}"
bundle="${P43_GATE2_BUNDLE:-p43-topology-gate2-v2.tar.gz}"
job_source="$study_root/gate2-source-$SLURM_JOB_ID"
mkdir "$job_source"
tar -xzf "$study_root/$bundle" -C "$job_source"
cd "$job_source"
if [[ "${P43_INSTALL_VERIFY_DEPS:-0}" == 1 ]]; then
  dependency_dir="$job_source/dependencies"
  mkdir "$dependency_dir"
  "$python_path" -m pip install --disable-pip-version-check --only-binary=:all: \
    --target "$dependency_dir" --retries 1 --timeout 12 --no-deps \
    'nlopt==2.9.1' 'shapely==2.1.2' 'triangle==20250106'
  if ! "$python_path" -c 'import pytest'; then
    "$python_path" -m pip install --disable-pip-version-check --only-binary=:all: \
      --target "$dependency_dir" --retries 1 --timeout 12 'pytest==8.4.2'
  fi
  export PYTHONPATH="$dependency_dir${PYTHONPATH:+:$PYTHONPATH}"
fi
"$python_path" -B -u -m project4_fem.project43_topology.gate2 "$study_root/gate2-results-$SLURM_JOB_ID" "$@"
