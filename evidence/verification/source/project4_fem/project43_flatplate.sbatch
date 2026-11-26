#!/bin/bash
#SBATCH -p htc
#SBATCH -J p43-flatplate
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 01:00:00
set -eu
cd "${SLURM_SUBMIT_DIR:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$HOME/.conda/envs/brinkman-rans/bin/python" -B -m project4_fem.project43_flatplate --output "flatplate-$SLURM_JOB_ID.json"
