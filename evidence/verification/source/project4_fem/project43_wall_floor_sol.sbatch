#!/bin/bash
#SBATCH --partition=htc
#SBATCH --job-name=p43-wall-floor
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=slurm-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
py="${P43_PYTHON:-$HOME/.conda/envs/brinkman-rans/bin/python}"
out="p43-wall-floor-${SLURM_JOB_ID:?}"
mkdir "$out"
sha256sum project4_fem/{fem,fillet,junction,project43_baseline,project43_brinkman,project43_rans_limit,project43_sst}.py project4_fem/next_rans_adjoint/sst.py project4_fem/project43_wall_floor_sol.sbatch > "$out/input-sha256.txt"
"$py" -B -c 'import sys,numpy,scipy; print(sys.version,numpy.__version__,scipy.__version__)' > "$out/environment.txt"
code=0
"$py" -B -u -m project4_fem.project43_sst --spacing-mm 0.5 --lump-destruction --refine-newton --pseudo-transient --wall-floor-fraction 0.5 --output "$out/report.json" > "$out/solve.log" 2>&1 || code=$?
if [[ -s "$out/report.json" ]]; then sha256sum "$out/report.json" > "$out/output-sha256.txt"; fi
exit "$code"
