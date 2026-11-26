#!/bin/bash -l
#SBATCH --job-name=p43-topology-recovery
#SBATCH --partition=htc
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:35:00
#SBATCH --array=0-2
#SBATCH --output=p43-topology-recovery-%A_%a.log

set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLBACKEND=Agg
study_root="${P43_STUDY_ROOT:-/scratch/$USER/p43-topology-study-20261125}"
dependency_dir="${P43_DEPS:-$study_root/gate2-source-2395294/dependencies}"
export PYTHONPATH="$dependency_dir${PYTHONPATH:+:$PYTHONPATH}"
python_path="${P43_PYTHON:-$HOME/.conda/envs/brinkman-rans/bin/python}"
job_source="$study_root/recovery-source-$SLURM_JOB_ID"
output="$study_root/recovery-results-$SLURM_ARRAY_JOB_ID"
mkdir "$job_source"
mkdir -p "$output"
tar -xzf "$SLURM_SUBMIT_DIR/p43-topology-recovery.tar.gz" -C "$job_source"
cd "$job_source"
anchor="$study_root/restoration-results-2395874/full"
case "$SLURM_ARRAY_TASK_ID" in
  0) /usr/bin/time -v "$python_path" -B -u -m project4_fem.project43_topology.reduced_authority \
    "$output/reduced" "$anchor/restoration_summary.json" --plot ;;
  1) /usr/bin/time -v "$python_path" -B -u -m project4_fem.project43_topology.full_recovery \
    "$anchor" "$output/slsqp" --algorithm slsqp ;;
  2) /usr/bin/time -v "$python_path" -B -u -m project4_fem.project43_topology.full_recovery \
    "$anchor" "$output/mma" --algorithm mma ;;
esac
