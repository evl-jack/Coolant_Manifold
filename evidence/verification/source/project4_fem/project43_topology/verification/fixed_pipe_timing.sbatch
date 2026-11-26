#!/bin/bash -l
#SBATCH --job-name=p43-sens-cost
#SBATCH --partition=htc
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00

set -euo pipefail
root="${1:?Pass the staging directory containing source}"
reference="${2:?Pass the archived stage-03 beta-8 original metrics.json}"
job="$SLURM_JOB_ID"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
export PYTHONPATH="$root/source:/scratch/$USER/p43-topology-study-20261125/gate2-source-2395294/dependencies"
python="${P43_PYTHON:-/scratch/$USER/p43-accelerated/cpu-py312-jax053/bin/python}"
cd "$root"

status=0
for case_name in pilot-128 scale-32 scale-128 scale-512; do
    workers="$SLURM_CPUS_PER_TASK"
    if [[ "$case_name" == pilot-128 ]]; then workers=1; fi
    printf 'Sol job %s, case %s, workers %d\n' "$job" "$case_name" "$workers"
    srun "$python" -B -u -m project4_fem.project43_topology.verification.fixed_pipe_timing \
        run "$case_name" "$root/timing-$job/$case_name" --reference-metrics "$reference" --workers "$workers" \
        || status=$?
done
if [[ "$status" -eq 0 ]]; then
    "$python" -B -u -m project4_fem.project43_topology.verification.fixed_pipe_timing \
        plot "$root/timing-$job" "$root/timing-$job/figures"
fi
exit "$status"
