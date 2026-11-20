# 3D manifold RANS on ASU Sol

This case uses the latest selected shape in `rans_solver/results/shape_best.csv`.
The main bore is revolved about x and the side branch about y. The 3.5 mm
main openings use 0.25 mm long lips. The main throat is 2.0 mm in diameter,
the side opening is 1.25 mm, and the branch extends into the main bore at the
junction. The STEP coordinates are in mm and the OpenFOAM mesh coordinates are
in meters. The 2D design area budget does not specify a 3D volume limit.

## VS Code

Open the project root and use the existing Python virtual environment.
Install the meshing dependency and regenerate the tagged surface if the 2D
shape changes. The surface contains the three port caps and the walls in one
closed STL, with coordinates in meters.

```powershell
.\.venv\Scripts\python.exe -m pip install -r sol_rans_3d/requirements.txt
.\.venv\Scripts\python.exe sol_rans_3d/build_mesh.py surface
```

OpenFOAM meshes this surface with `blockMesh`, `surfaceFeatureExtract`, and
`snappyHexMesh`. The coarse and fine background grids use uniform nominal
cell sizes of 0.09 and 0.06 mm. Uniform sizing avoids nonplanar faces at local
refinement transitions, which the full mesh check flagged. The snap step
uses the geometry's existing lips, branch junction, and outlet caps. The
legacy tetrahedral meshes can still be generated with `coarse` or `fine` as
the argument to `build_mesh.py`.

## ASU Sol

Connect to the ASU VPN and use [ASU's SSH instructions](https://docs.rc.asu.edu/new-user-guide/)
or [Open OnDemand](https://docs.rc.asu.edu/web-portal/). Upload the bundle to your
Sol scratch directory and unpack it. Use your assigned Slurm account or QoS if
your project or class requires one.

```bash
cd "/scratch/$USER"
tar -xzf sol_rans_3d.tar.gz
cd sol_rans_3d
container_job=$(sbatch --parsable pull_openfoam.sbatch)
coarse_job=$(sbatch --parsable --dependency=afterok:"${container_job%%;*}" submit.sbatch coarse)
squeue -u "$USER"
```

This uses the [OpenCFD v2312 image](https://hub.docker.com/r/opencfd/openfoam-default/tags)
with [Apptainer](https://docs.rc.asu.edu/apptainer/). If a compatible OpenCFD
module is available, submit directly using
`sbatch --export=ALL,FOAM_MODULE=<module> submit.sbatch coarse`.
The solver job defaults to one node with 16 MPI ranks.

Use the tagged STL on native OpenFOAM v2512 with `MESH_MODE=hex`. Upload the
updated case directory, scripts, and `meshes/manifold.stl` before submitting.
For a mesh-only check, add `MESH_CHECK_ONLY=1` and request one task. Use
`fine` as the final argument and `--mem=20G` to check the finer grid.

```bash
sbatch --ntasks=1 --mem=12G --time=00:12:00 --export=ALL,FOAM_MODULE=openfoam/2512,MESH_MODE=hex,MESH_CHECK_ONLY=1 submit.sbatch coarse
```

If concave cells are the only reported issue, copy the run's `constant/polyMesh`
directory to Windows and use the project Python environment to smooth the
flagged cells. The least-squares fit moves nearby interior vertices and
allows small wall adjustments. Port vertices stay fixed. The script writes
a replacement `points` file and reports the maximum wall movement.

```powershell
.\.venv\Scripts\python.exe sol_rans_3d/repair_hex_mesh.py .\polyMesh 0.09 .\points.repaired
```

On Sol, keep the initial mesh for comparison. Copy its `polyMesh` to a
separate directory, place `points.repaired` in that directory as `points`,
and use the repaired directory in `MESH_PREBUILT`. This avoids regenerating
the grid and runs the same `checkMesh -allGeometry -allTopology` gate before
starting `simpleFoam`.

```bash
sbatch --export=ALL,FOAM_MODULE=openfoam/2512,MESH_PREBUILT="$PWD/repaired-coarse/constant/polyMesh" submit.sbatch coarse
```

Review `runs/coarse-<job ID>/log.checkMesh`, `log.simpleFoam`, and `summary.json`
before solving on the fine mesh. The job reconstructs the latest field
for ParaView. Open `manifold.foam` in that run directory.
The native module launch follows the Sol OpenMPI example. [ASU MPI launchers](https://docs.rc.asu.edu/mpi/)

### Staged upwind convergence test

Add `STAGED_UPWIND=1` to the Slurm export list to start a case with first-order
upwind convection for `U`, `k`, and `omega`. After 300 iterations the run
compares the first 10 and last 20 median initial residuals for `p`, `U`, `k`,
and `omega`. Velocity uses the largest component residual at each iteration.
The run restarts from its saved fields with linear upwind for all three
transport terms if every residual falls by at least 100 times. The final
iteration cap and convergence tolerances come from the case dictionaries.
`upwind_switch.json` records the decision and both median residuals. If the
case switches, `log.upwind` holds the first stage and `log.simpleFoam` holds
the final stage.

### Transient 3D simulation

Use `submit_urans.sbatch` with a completed steady-case directory and a saved,
reconstructed field time. The batch job copies those fields to a new physical
time zero and uses the repaired mesh, inlet and outlet conditions, and SST
model from the seed case.

```bash
sbatch submit_urans.sbatch path/to/steady-case 300 0.02 0.01
```

The last two arguments give the end time and averaging start in seconds.
The job runs `pimpleFoam` with backward time differencing and adaptive
time steps targeting a maximum Courant number of 1. Momentum uses linear
upwind. Turbulence transport uses the seed’s upwind discretization. The
job writes a flow field every 2 ms, saves `UMean` and `pMean` starting at
the averaging time, and reconstructs the saved fields for ParaView.

`history.csv` records all monitored time steps. `summary.json` compares
time-weighted flow, pressure drop and reverse flow over the average window
and its two halves. The pressure and flow ranges help assess whether the
selected window is long enough to represent the transient behavior.

Open `manifold.foam` in a reconstructed run to play the time series.
`render_urans.py` renders the saved times with a 3D view of the walls,
velocity on the central slice, and instantaneous streamlines.

```bash
pvbatch render_urans.py path/to/run/manifold.foam frames
```

## What to compare

The [steady incompressible case](https://doc.openfoam.com/2312/tools/processing/solvers/rtm/incompressible/simpleFoam/)
uses water with density 998.2 kg/m³, an inlet speed of 2 m/s, two outlets at zero
gauge pressure, and explicit no-slip walls with k–omega SST. Compare the
coarse and fine outlet flows, static pressure drop, and wall y+ before
interpreting the results. The target is `Q_side / Q_main = 0.75`. The
outlet pressures determine the actual split. The tetrahedral meshes use wall
refinement. Inspect y+ and add a wall layer mesh if near-wall resolution is
insufficient for the quantity being assessed.
