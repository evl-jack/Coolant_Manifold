#!/bin/bash
#SBATCH --partition=htc
#SBATCH --job-name=p43-v45
#SBATCH --array=0-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
set -eu
cd "${SLURM_SUBMIT_DIR:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
py="$HOME/.conda/envs/brinkman-rans/bin/python"
out="verification-${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}"
mkdir "$out"
case "$SLURM_ARRAY_TASK_ID" in
  0) "$py" -B -m project4_fem.project43_leakage "/scratch/$USER/brinkman-run-2349118" > "$out/leakage.json" ;;
  1|2) floors=(0 0 0.5)
       "$py" -B -u -m project4_fem.project43_sst --spacing-mm 0.5 --lump-destruction --refine-newton --pseudo-transient --wall-floor-fraction "${floors[$SLURM_ARRAY_TASK_ID]}" --output "$out/report.json" > "$out/solve.log" 2>&1 ;;
esac
