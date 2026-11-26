#!/bin/bash
set -euo pipefail
export PYTHONPATH="$PWD/source" MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1
py="/scratch/$USER/p43-accelerated/cpu-py312-jax053/bin/python"
"$py" -B -m project4_fem.project43_topology.verification.plot_fixed_pipe_adjoint adjoint-coordinate-supplement.json adjoint-coordinate-evidence
sacct -j "$(cut -d= -f2 jobs.txt | paste -sd,)" --parsable2 --format=JobID,JobName,State,ExitCode,NodeList,Elapsed,TotalCPU,MaxRSS,ReqMem > slurm_accounting.txt
zip -q -r ../Fixed_Pipe_Sol_Evidence_stage.zip source input mesh-results mesh-evidence adjoint-results adjoint-evidence adjoint-coordinate-evidence kkt-results timing-2408534 *.log jobs.txt slurm_accounting.txt launch_sol_audit.sh coordinate_supplement.py adjoint-coordinate-supplement.json bundle_evidence.sh -x "*/__pycache__/*"
mv ../Fixed_Pipe_Sol_Evidence_stage.zip ../Fixed_Pipe_Sol_Evidence.zip
unzip -tq ../Fixed_Pipe_Sol_Evidence.zip
sha256sum ../Fixed_Pipe_Sol_Evidence.zip
