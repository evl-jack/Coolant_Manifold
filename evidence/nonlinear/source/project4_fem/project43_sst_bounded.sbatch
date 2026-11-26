#!/bin/bash
#SBATCH --partition=htc
#SBATCH --job-name=p43-sst-bounded
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out
set -eu
cd "${SLURM_SUBMIT_DIR:?}"
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
out="run-${SLURM_JOB_ID:?}"
mkdir "$out"
sha256sum project4_fem/*.py project4_fem/next_rans_adjoint/sst.py > "$out/source.sha256"
"$HOME/.conda/envs/brinkman-rans/bin/python" -B -u -m project4_fem.project43_sst_bounded \
    "${1:?Supply the sweep-200 checkpoint}" "$out"
