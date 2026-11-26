#!/bin/bash
#SBATCH --partition=htc
#SBATCH --job-name=p43-verify
#SBATCH --array=0-2
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
set -eu
cd "${SLURM_SUBMIT_DIR:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cases=(geometry laminar brinkman)
case_name=${cases[$SLURM_ARRAY_TASK_ID]}
"$HOME/.conda/envs/brinkman-rans/bin/python" -B -m project4_fem.project43_verify_three "$case_name" > "verify-${SLURM_ARRAY_JOB_ID}-$case_name.json"
