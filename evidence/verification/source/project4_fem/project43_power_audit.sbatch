#!/bin/bash
#SBATCH --partition=htc
#SBATCH --job-name=p43-power-audit
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
set -eu
cd "${SLURM_SUBMIT_DIR:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$HOME/.conda/envs/brinkman-rans/bin/python" -B -m project4_fem.project43_power_audit "${1:?split feasible directory}" "power-audit-$SLURM_JOB_ID.json"
