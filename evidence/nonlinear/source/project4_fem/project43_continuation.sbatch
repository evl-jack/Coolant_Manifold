#!/bin/bash -l
#SBATCH -p htc
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 01:15:00
#SBATCH -J p43-sst-cont
set -eu
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd "${SLURM_SUBMIT_DIR:?}"
checkpoint="${1:?Provide the archived beta=1 sweep-200 checkpoint}"
sha256sum "$checkpoint" project4_fem/project43_sst_*.py > "sst-cont-$SLURM_JOB_ID.sha256"
"$HOME/.conda/envs/brinkman-rans/bin/python" -B -m project4_fem.project43_sst_continuation \
    "$checkpoint" "sst-cont-$SLURM_JOB_ID"
