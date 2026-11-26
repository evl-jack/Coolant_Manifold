#!/bin/bash -l
#SBATCH -p htc
#SBATCH -J p43-laminar-kkt
#SBATCH -t 00:10:00
#SBATCH -c 1
#SBATCH --mem=4G

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the extracted audit directory}"
archive="${1:?Pass the absolute path of lam-pair-01}"
out="$PWD/results-$SLURM_JOB_ID"
mkdir "$out"
sha256sum -c package.sha256 > "$out/package-integrity.log"
cp package.sha256 "$out/audit-source.sha256"
(cd "$archive"; sha256sum pair.json project4_fem/{lam_pair,project43_baseline,project43_brinkman,project43_power,fem,junction,fillet,project43_explicit}.py) > "$out/archived-source.sha256"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONPATH="$archive"
python="$HOME/.conda/envs/brinkman-rans/bin/python"
"$python" -c 'import numpy, scipy, sys; print(sys.version, numpy.__version__, scipy.__version__)' > "$out/environment.txt"
status=0
"$python" -B project4_fem/project43_laminar_endpoint_audit.py "$archive" > "$out/audit.json" 2> "$out/audit.log" || status=$?
printf '%s\n' "$status" > "$out/exit_code.txt"
exit "$status"
