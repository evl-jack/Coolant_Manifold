# Project 4.3 SST handoff · 24 November 2026

**A new collaborator can start with the saved flat-plate state.** We have a working finite-element SST implementation and a converged low-Reynolds-number junction case. At the turbulent flat-plate operating point, two initial guesses produced substantially different skin friction and both stalled. The immediate engineering question is whether a different forward iteration can converge the *same* discrete equations well enough to predict wall shear, junction flow split, and pumping power.

## Physical setup and code

The project is a 2D inlet–main outlet–side outlet manifold. The design objective is minimum port pumping power per unit depth at prescribed inlet flow, subject to a side share between 49.9% and 50.1% and a geometric area limit of 83 mm². The Re = 50 laminar design and its successful shoulder derivative check (Sol 2350814) provide a workflow reference. The turbulent design requires its own operating point, inlet turbulence conditions, baseline, and derivative check.

| File | Start here for |
|---|---|
| [next_rans_adjoint/sst.py](next_rans_adjoint/sst.py) | SST-2003m closure, eddy viscosity, production limiter, blending, conservative scalar transport, streamline stabilization, complex-step element Jacobian |
| [project43_flatplate.py](project43_flatplate.py) | Flat-plate mesh, essential conditions, initial state, lumped destruction, and wall gradient |
| [project43_sst.py](project43_sst.py) | Pseudo-transient iteration, turbulence positivity damping, residual term decomposition, and frozen-geometry Brinkman junction |
| [fem.py](fem.py) | Taylor–Hood mesh, quadratic velocity interpolation, linear pressure, geometry and boundary-edge ordering |
| [PROJECT43_WALL_FLOOR.md](PROJECT43_WALL_FLOOR.md) | Matched junction runs and sensitivity to regularizing the geometric wall distance |

The fluid is steady, incompressible and constant-density. The FEM unknowns are quadratic velocity, linear modified pressure, and linear turbulence scalars. The external-flow test applies no Brinkman force. The plate is fully turbulent in the modeling sense: SST contains no transition model here.

| Plate input | Implementation and interpretation |
|---|---|
| Domain and mesh | Inlet at x = −1/3, leading edge x = 0, plate end x = 2, top y = 1. Re = 5 × 10⁶ per unit length, U∞ = 1, 3,956 triangles, first scalar-node height 4 × 10⁻⁶. |
| Inlet and top | Prescribed (u,v) = (1,0), k = 1.125/Re = 2.25 × 10⁻⁷, ω = 125. The turbulence values correspond to the benchmark at this Reynolds number. |
| Lower boundary | Upstream symmetry imposes zero normal velocity. On the plate impose no slip, k = 0 and ω = 60ν/(0.075 Δy₁²), equal to 10⁷ in these units. Only the plate supplies the wall distance for SST. Positive k is required at free interior nodes; zero is correct at the plate. |
| Outlet | Natural total traction for momentum, outgoing turbulence flux with interior trace and natural diffusive flux. Check for reverse flow before interpreting this boundary. |
| Observable | Cf = 2(∂u/∂y + ∂v/∂x)/Re at x = 0.97008. The recorded comparator is 0.00268299 from the CFL3D 273 × 193 SST-Vm result. |
| Differences to reference | TMR uses compressible Mach 0.2 SST-Vm with different outer boundary treatments and grid family. Our incompressible SST-2003m screen measures aggregate discrepancy. |

## Evidence from Sol

Paths below are relative to `/scratch/jtfoste6/`. A trial Cf printed in stdout is diagnostic. The pilot driver records a null accepted Cf and exit code 2 on solver failure.

| Job and working directory | Finding |
|---|---|
| 2349247, tasks 1 and 2 (junction) | At Re = 50 and 0.5 mm spacing, the zero-distance-floor case reached residual 2.61 × 10⁻¹³, side share 4.57918%, and static pressure loss 0.164719 Pa. Applying a 0.5-local-altitude distance floor yielded 1.64 × 10⁻¹¹ against the original 10⁻¹¹ threshold. This is a different geometry and regime from the turbulent plate. |
| 2349550, `flatplate-pilot/` | Seventh-power-root velocity start with prescribed boundary-layer turbulence profiles. The solver exhausted its pseudo-time budget at residual 0.55996. |
| 2350393, `sst-diag-ZOeWNg/` | The free ω degree of freedom at node 870, (0.603738, 0.00061227), reached ≈3.56 × 10⁻¹³ while its proposed change was about −18,000 to −28,000. The global positivity factor was of order 10⁻¹⁷. |
| 2350508, `sst-row-BAsXWh/` | Saved `precollapse.npz` when node 870 first limited the update. Its ω was 3.69 versus 10.67–60.60 at neighboring nodes. Local directional derivatives were 3.957311405 (Jacobian) and 3.957311543 / 3.957311367 (centered differences). The largest local residual contributions were diffusion −0.07203 and production −0.02141. This checks one row and direction; the full coupled update depends on other rows. |
| Same checkpoint | Independent six-node quadratic reconstruction reproduced Cf = 0.06101222434. The original initial field gave 0.07979; the final iterate gave 0.0610061. The large shear is present in the velocity field. |
| 2350713, `sst-wallstart-v1NisD/` | A linear viscous near-wall velocity start used local Cf ≈ 0.0592 Reₓ⁻¹ᐟ⁵. The last Cf was 0.00281075 (4.76% high), residual 0.97281, with k at node 805 near the leading edge limiting the update. The correlation already initializes Cf close to the comparator. Solver exit 2. |
| 2351110, `sst-independent-01/` | Submitted with uniform positive free k and ω, the prescribed wall traces, and downstream u = tanh(y/0.02) before applying boundary values. The wall slope is independent of the benchmark Cf. See stdout and JSON for the terminal result; the job ID alone records submission. |

Read the job-specific JSON and `slurm-<jobid>.out` together. The staged project folders contain the diagnostic variants used for each job. The checkpoint stores `state`, `free`, and `update`; use the mesh constructor staged with that job to decode node numbers. The reported maximum residual uses separate initial block scales with a lower bound of one. Across different starts, compare unscaled block residuals and physical outputs as well as normalized values.

The following Sol commands retrieve the existing runs and checkpoint diagnostics.

```bash
sacct -X -j 2349550,2350393,2350508,2350713,2351110 -o JobID,State,Elapsed,ExitCode
cd /scratch/jtfoste6/sst-row-BAsXWh
grep -E 'row tangent|omega balance|neighbors' slurm-2350508.out
ls -lh precollapse.npz
cd /scratch/jtfoste6/sst-independent-01
tail -12 slurm-2351110.out
tail -15 flatplate-2351110.json
```

## What to work on independently

1. Trace the leading-edge k balance in the correlation-start run. Inspect the free interior node immediately over the symmetry/no-slip junction, its neighbors, and the local advection, diffusion and source contributions. Determine whether the update creates an undershoot or amplifies one already present.
2. Test a segregated forward update. Update velocity/pressure and scalar fields in separate blocks, lag production and treat destruction implicitly. Recompute the complete original steady residual after each outer sweep. At convergence, retain every turbulence coupling in the physical Jacobian. If the converged residual itself changes through stabilization or limiting, differentiate that actual residual.
3. Compare the three downstream shear stations chosen **before** evaluation, outer velocity and pressure gradient against the TMR case. The stock pilot logs one shear station; adding the other two is required for a three-station claim. Trial engineering screen proposed in discussion: all three local Cf values within 5% of reference, less than 1% drift over ten accepted updates, materially reduced block residuals, and admissible free k and ω. Record actual values and acceptance history. These are project criteria, not published tolerances.
4. On the fixed turbulent baffle, check power, split, mass conservation, leakage and explicit-wall comparison, including dependence on independent initial states. For optimization, perturb one shoulder direction twice at decreasing step sizes; reconverge the flow and recompute wall distance for each perturbation. Compare signs and magnitudes of both objective and flow-split derivatives. A split constraint margin smaller than the measured discretization or representation error calls for an independent final geometry check.

## What “frozen” means for Brinkman

The current reference fixes a binary, mesh-aligned interface, includes finite momentum drag, removes scalar residuals in solid cells and sets wall traces on known nodes. During an evolving-density solve these switches would introduce discontinuities. The proposed extension uses γ = 1 in fluid, 0 in solid, and s = (1−γ)². Add αmax s u/ρ in the momentum residual, λk s k in the k residual and λω s(ω−ωtarget) in the ω residual, with consistent nondimensionalization. These scalar penalties, their strengths, and a design-dependent distance treatment are modeling proposals requiring verification. Lagging coefficients inside a forward iteration can retain the same steady solution; freezing turbulence in the final design derivative drops dependencies. A wall-distance floor regularizes the blending length. A scalar floor alters the treatment of a different quantity and would require its own sensitivity and activation diagnostics.

## Source material

- Published flat-plate conditions and CFL3D/FUN3D results, including first-order upwinding of turbulence advection and the grid-dependence of skin friction. [TMR case](https://tmbwg.github.io/turbmodels/flatplate.html) · [SST results and data](https://tmbwg.github.io/turbmodels/flatplate_sst.html) · [SST definitions and wall conditions](https://tmbwg.github.io/turbmodels/sst.html)
- Carrier et al. used finite-element SST with automatic wall treatment and assessed engineering-output accuracy on the same plate. Their approach is a reference for wall treatment and mesh design. [International Journal for Numerical Methods in Fluids](https://doi.org/10.1002/fld.4905)
- Dilgen et al. used related finite-volume turbulent topology optimization, explicit production/implicit dissipation handling, and discrete adjoints. Their k–ω treatment informs an iteration strategy; our SST FEM implementation needs its own checks. [Computer Methods in Applied Mechanics and Engineering](https://doi.org/10.1016/j.cma.2017.11.029)

## Sol implementation update

Independent-start control 2351110 showed Cf 0.00002035 and scaled residual 0.99617 at iteration 68. This is an intermediate observation.

Segregated pilot 2351184 ran in /scratch/jtfoste6/sst-segregated-01. The exact executable source is project4_fem/segregated.py there. It imports the instrumented independent-start flat-plate driver from its own directory. Launch from that directory using the brinkman-rans Python environment and module project4_fem.segregated with a fresh JSON output filename. Each sweep updates the mixed velocity/pressure block, then k and omega with lagged transport coefficients, explicit production and implicit lumped destruction. Each scalar has its own positivity limiter. The original coupled residual is evaluated after each sweep. This pilot retains the original scalar transport discretization. Bounded advection, moving density and a turbulent optimization adjoint have not been implemented in this pilot. Monitor slurm-2351184.out for iteration, residual, Cf, and the separate scalar step factors. The first two iterations executed. Acceptance is pending. JSON includes an explicit convergence flag; check it even when Slurm reports success.

## Scalar transport experiment

Job 2351184 reached sweep 99 with original scaled residual 0.99146 and positivity factors below 1e-91 for both scalar blocks. Separating the updates alone did not cure the stall. Job 2351259 in /scratch/jtfoste6/sst-bounded-01 adds conservative pairwise graph diffusion to the scalar residual and its lagged scalar operator. The stored residual metric is the original unstabilized residual. It cannot certify convergence of the modified equations. Inspect the scalar update factors and friction history; report the stabilized residual and stabilization contribution separately before using this variant for optimization. There is no hard turbulence floor.
