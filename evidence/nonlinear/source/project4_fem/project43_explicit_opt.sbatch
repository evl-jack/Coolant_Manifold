#!/bin/bash
#SBATCH --partition=htc
#SBATCH --job-name=p43-wall-opt
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
set -eu
cd "${SLURM_SUBMIT_DIR:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$HOME/.conda/envs/brinkman-rans/bin/python" -B -u -m project4_fem.project43_explicit_opt "wall-opt-$SLURM_JOB_ID"
