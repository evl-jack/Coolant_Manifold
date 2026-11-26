#!/bin/bash -l
#SBATCH --job-name=p43-hybrid-diagnostic
#SBATCH --partition=htc
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=01:30:00
#SBATCH --output=p43-hybrid-%j.log

# From the directory containing this file and hybrid-diagnostics.tar.gz:
# verify=$(sbatch --parsable hybrid_diagnostics.sbatch)
# seed=$(sbatch --parsable --dependency=afterok:$verify --export=ALL,P43_HYBRID_MODE=seed,P43_HYBRID_VERIFY_JOB=$verify hybrid_diagnostics.sbatch)
# For each of model-a and model-b, submit with --dependency=afterok:$seed and
# --export=ALL,P43_HYBRID_MODE=<model>,P43_HYBRID_VERIFY_JOB=$verify,P43_HYBRID_SEED_JOB=$seed
# Supply the same P43_HYBRID_MAX_STAGES to both model jobs to cap the matched schedule.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
root="/scratch/$USER/p43-topology-study-20261125"
python_path="$HOME/.conda/envs/brinkman-rans/bin/python"
bundle="$SLURM_SUBMIT_DIR/hybrid-diagnostics.tar.gz"
source_dir="$root/hybrid-source-$SLURM_JOB_ID"
mode="${P43_HYBRID_MODE:-verify}"
results="$root/hybrid-$mode-$SLURM_JOB_ID"
export PYTHONPATH="$root/gate2-source-2395294/dependencies"

mkdir "$source_dir" "$results"
tar -xzf "$bundle" -C "$source_dir"
sha256sum "$bundle" | cut -d ' ' -f 1 > "$results/source.sha256"
cd "$source_dir"

run() {
    module="$1"
    shift
    /usr/bin/time -v "$python_path" -B -u -m "project4_fem.project43_topology.$module" "$@"
}

if [[ "$mode" == verify ]]; then
    run verify_hybrid "$results/hybrid_verification.json"
    exit
fi

verify_dir="$root/hybrid-verify-${P43_HYBRID_VERIFY_JOB:?Set the passed verification job ID}"
cmp "$verify_dir/source.sha256" "$results/source.sha256"
"$python_path" -B -c 'import json, sys; sys.exit(0 if json.load(open(sys.argv[1]))["passed"] else 1)' \
    "$verify_dir/hybrid_verification.json"

arguments=(--mode "$mode" --route coupled --mma-evaluations 150 --restoration-evaluations 250)
if [[ -n "${P43_HYBRID_MAX_STAGES:-}" ]]; then
    arguments+=(--max-stages "$P43_HYBRID_MAX_STAGES")
fi

case "$mode" in
    seed) ;;
    model-a|model-b)
        seed_dir="$root/hybrid-seed-${P43_HYBRID_SEED_JOB:?Set the common seed job ID}"
        cmp "$seed_dir/source.sha256" "$results/source.sha256"
        arguments+=(--shared-seed "$seed_dir/diagnostic/seed.json") ;;
    *) echo "Unsupported mode: $mode" >&2; exit 64 ;;
esac

run hybrid_driver "$results/diagnostic" "${arguments[@]}"
