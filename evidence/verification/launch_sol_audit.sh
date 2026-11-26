#!/bin/bash
set -euo pipefail
root="$(cd "$(dirname "$0")" && pwd)"
source_dir="$root/source"
audit="$source_dir/project4_fem/project43_topology/verification"
python="/scratch/$USER/p43-accelerated/cpu-py312-jax053/bin/python"
deps="/scratch/$USER/p43-topology-study-20261125/gate2-source-2395294/dependencies"
export PYTHONPATH="$source_dir:$deps"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
"$python" -c 'import scipy, numpy, matplotlib, shapely, triangle; print("Sol dependencies available")'
sha256sum "$root/input/stage_03_beta_8/"{original,selected}/fields.npz "$root/input/smoothed_geometry.json"
find "$source_dir" -type f -name '*.sbatch' -exec sed -i 's/\r$//' {} +
cd "$root"
mesh=$(sbatch --parsable -o "$root/mesh-%A-%a.log" "$audit/mesh_cases.sbatch" "$root" "$root/mesh-results")
echo "mesh=$mesh" | tee -a jobs.txt
adj=$(sbatch --parsable -o "$root/adjoint-%A-%a.log" "$audit/fixed_pipe_adjoint_audit.sbatch" "$source_dir" "$root/input" "$root/adjoint-results")
echo "adjoint=$adj" | tee -a jobs.txt
kkt=$(sbatch --parsable -J p43-kkt-audit -p htc -c 1 --mem=6G -t 00:40:00 -o "$root/kkt-%j.log" --wrap "'$python' -B -u -m project4_fem.project43_topology.verification.fixed_pipe_kkt '$root/input' '$root/kkt-results' --replay-history")
echo "kkt=$kkt" | tee -a jobs.txt
cost=$(sbatch --parsable -o "$root/timing-%j.log" "$audit/fixed_pipe_timing.sbatch" "$root" "$root/input/stage_03_beta_8/original/metrics.json")
echo "timing=$cost" | tee -a jobs.txt
report=$(sbatch --parsable -J p43-mesh-report -p htc -c 1 --mem=8G -t 00:20:00 -d "afterany:${mesh%%;*}" -o "$root/report-%j.log" --wrap "'$python' -B -u -m project4_fem.project43_topology.verification.mesh_report '$root/mesh-results' '$root/mesh-evidence'")
echo "mesh_report=$report" | tee -a jobs.txt
