#!/bin/bash
#SBATCH --partition=htc
#SBATCH --job-name=p43-sst-check
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=slurm-%j.out
set -eu
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "${SLURM_SUBMIT_DIR:?}"
py="$HOME/.conda/envs/brinkman-rans/bin/python"
test -x "$py"
out="run-${SLURM_JOB_ID:?}"
mkdir "$out"
sha256sum project4_fem/{fem,fillet,junction,project43_baseline,project43_brinkman,project43_rans_limit,project43_sst}.py project4_fem/next_rans_adjoint/sst.py project4_fem/project43_sst_sol.sbatch > "$out/input-sha256.txt"
"$py" -B -c 'import sys,numpy,scipy; print(sys.version,numpy.__version__,scipy.__version__)' > "$out/environment.txt"
case_result=0
for variant in control lumped; do
    options=()
    if [[ $variant == lumped ]]; then options+=(--lump-destruction --refine-newton); fi
    code=0
    "$py" -B -u -m project4_fem.project43_sst --spacing-mm 0.5 "${options[@]}" --output "$out/$variant.json" > "$out/$variant.log" 2>&1 || code=$?
    printf '%s %d\n' "$variant" "$code" >> "$out/exit-codes.txt"
    [[ $code == 0 || ( $code == 2 && -s "$out/$variant.json" ) ]] || exit "$code"
    if [[ $variant == lumped ]]; then case_result=$code; fi
done
exit "$case_result"
