#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
python="${P43_PYTHON:-python3}"
export PYTHONPATH="$PWD/source${P43_DEPENDENCIES:+:$P43_DEPENDENCIES}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
module=project4_fem.project43_topology.verification
out="$(mktemp -d "$PWD/rebuilt-XXXXXX")"
"$python" -B -m "$module.plot_fixed_pipe_adjoint" adjoint-evidence/adjoint_merged.json "$out/adjoint"
"$python" -B -m "$module.plot_fixed_pipe_adjoint" adjoint-coordinate-supplement.json "$out/adjoint-coordinate"
"$python" -B -m "$module.mesh_report" mesh-results "$out/mesh"
"$python" -B -m "$module.fixed_pipe_timing" plot timing-2408534 "$out/timing"
"$python" -B -u -m "$module.fixed_pipe_kkt" input "$out/kkt" --replay-history
printf '%s\n' "$out"
