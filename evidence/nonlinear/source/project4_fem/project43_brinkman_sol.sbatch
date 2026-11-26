#!/bin/bash
#SBATCH --partition=htc
#SBATCH --job-name=p43-lam-brinkman
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=slurm-%j.out

# Submit from the directory containing project4_fem.
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
py="${P43_PYTHON:-$HOME/.conda/envs/brinkman-rans/bin/python}"
test -x "$py"

run_dir="project43-brinkman-${SLURM_JOB_ID:?}"
mkdir "$run_dir"
sha256sum project4_fem/{fem,fillet,junction,project43_baseline,project43_brinkman}.py \
    project4_fem/project43_brinkman_sol.sbatch > "$run_dir/input-sha256.txt"
"$py" -B -c 'import platform,sys,numpy,scipy; print(platform.platform(), sys.version, numpy.__version__, scipy.__version__)' \
    > "$run_dir/environment.txt"

# The nominal solve includes two centered re-solves to verify the log-penalty adjoint.
"$py" -B -u -m project4_fem.project43_brinkman --spacing-mm 0.25 --reynolds 50 \
    --alpha-max 1000000 --output "$run_dir/nominal" > "$run_dir/solve.log" 2>&1
sha256sum "$run_dir/nominal/report.json" "$run_dir/nominal/field.npz" > "$run_dir/output-sha256.txt"
