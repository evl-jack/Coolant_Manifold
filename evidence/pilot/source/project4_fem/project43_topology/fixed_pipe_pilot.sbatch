#!/bin/bash -l
#SBATCH --job-name=p43-fixed-pipe
#SBATCH --partition=htc
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00

set -euo pipefail
root="${1:?Pass the directory containing the staged source and launch command}"
job="$SLURM_JOB_ID"
output="$root/run-$job"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
export PYTHONPATH="$root/source:/scratch/$USER/p43-topology-study-20261125/gate2-source-2395294/dependencies"
python="/scratch/$USER/p43-accelerated/cpu-py312-jax053/bin/python"
cd "$root"

status=0
"$python" -B -u -m project4_fem.project43_topology.fixed_pipe_pilot "$output" || status=$?
printf 'Sol job %s, Python exit status %s\n' "$job" "$status"
tar --exclude='__pycache__' -czf "$root/Fixed_Pipe_Pilot_${job}.tar.gz.part" \
    source "run-$job" "sol-$job.log" launch.txt
mv "$root/Fixed_Pipe_Pilot_${job}.tar.gz.part" "$root/Fixed_Pipe_Pilot_${job}.tar.gz"
exit "$status"
