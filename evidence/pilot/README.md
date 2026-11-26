# Fixed-pipe Brinkman transition pilot

The summaries compare each selected immersed result with the original raw density replayed on the same mesh at the same q and beta. Threshold contours are diagnostic geometry. Effective hydraulic resistance is mechanical port power per unit depth divided by squared inlet flow per unit depth. Wall turns are sampled over 0.25 mm secants near the transition.

The five Windows benchmark evaluations replay Sol MMA unique evaluations 2 to 6 at beta 0, after warming up on the original design. Both runs used one BLAS thread. Sol times are elapsed differences between flushed callbacks. The local benchmark records peak Python working set.

## Local launch

```powershell
$env:OPENBLAS_NUM_THREADS='1'; $env:OMP_NUM_THREADS='1'; $env:MKL_NUM_THREADS='1'; $env:MPLBACKEND='Agg'
& .venv/Scripts/python.exe -B -u -m project4_fem.project43_topology.fixed_pipe_pilot --local output/fixed_pipe_pilot_local
```

## Earlier Sol launch, job 2408164

```bash
sbatch --parsable --output=/scratch/jtfoste6/p43-fixed-pipe-pilot-v1/sol-%j.log /scratch/jtfoste6/p43-fixed-pipe-pilot-v1/source/project4_fem/project43_topology/fixed_pipe_pilot.sbatch /scratch/jtfoste6/p43-fixed-pipe-pilot-v1
```

The Sol source, log and data archive is `/scratch/jtfoste6/p43-fixed-pipe-pilot-v1/Fixed_Pipe_Pilot_2408164.tar.gz`.
