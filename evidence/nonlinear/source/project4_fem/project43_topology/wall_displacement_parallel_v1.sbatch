#!/bin/bash -l
#SBATCH --partition=htc
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=wall-displacement-v1-%j.log
set -eu
root="$P43_WALL_ROOT"
cd "$root/source-linear-v1"
export PYTHONPATH="$PWD:/scratch/$USER/p43-topology-study-20261125/gate2-source-2395294/dependencies"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
python="/scratch/$USER/p43-accelerated/cpu-py312-jax053/bin/python"
package=project4_fem.project43_topology
output="$root/$P43_TAG-$SLURM_JOB_ID"
case "$P43_MODE" in
    geometry)
        exec "$python" -B -u -m "$package.wall_displacement_archive_audit_v1" "$P43_LANE" "$output.json";;
    evidence)
        exec "$python" -B -u -m "$package.wall_displacement_evidence_audit_v1" "$P43_LANE" --live-snapshot > "$output.json";;
    mesh)
        exec "$python" -B -u -m "$package.wall_displacement_physical_v1" "$root/$P43_TAG-$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID" \
            --optimization "$P43_LANE" --mesh-level "$SLURM_ARRAY_TASK_ID";;
    physical)
        exec "$python" -B -u -m "$package.wall_displacement_physical_v1" "$output" --optimization "$P43_LANE" --level-records \
            "$root/$P43_MESH_TAG-$P43_MESH_JOB-0/fitted_mesh_level.json" \
            "$root/$P43_MESH_TAG-$P43_MESH_JOB-1/fitted_mesh_level.json" \
            "$root/$P43_MESH_TAG-$P43_MESH_JOB-2/fitted_mesh_level.json";;
    matched)
        jq -e .results.physical_screen_passed "$P43_PHYSICAL_REPORT" >/dev/null || exit 2
        exec "$python" -B -u -m "$package.wall_displacement_conventional_v1" "$output" --candidate "$P43_PHYSICAL_REPORT";;
    export)
        exec "$python" -B -u -m "$package.wall_displacement_export_v1" "$output" --lane "$P43_LANE" --live-snapshot;;
esac
