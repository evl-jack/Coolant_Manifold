# Fixed-pipe Brinkman pilot verification

The paired mesh screen verified a pumping-power saving for the raw reconstructed wall while maintaining the original fitted side-flow fraction.

## Discrete adjoint

For a converged discrete state satisfying $R(z,m)=0$, solve $R_z^T\lambda_F=F_z$ and evaluate $\nabla_mF=F_m-R_m^T\lambda_F$. The fixed-port responses have zero direct design term. The transpose chain passes through the cell-to-quadrature map, rational Brinkman interpolation, projection and anchored Helmholtz filter. The nonlinear tangent differentiates both velocity factors in convection. See [the derivation](source/project4_fem/project43_topology/verification/ADJOINT_METHOD.md).

Finite differences passed for 15 of 15 tested directions. The initial Taylor order screen passed for 9 directions. The larger-step follow-up passed for 6 of 6 previously flagged directions and retained the original Taylor floor. Measured follow-up orders ran from 2.0010 to 2.0060.

The full 128-variable centered gradient used 257 forward solves. Its normalized-power and side-fraction infinity errors were 2.376e-12 and 1.099e-12.

| Cost case | Controls | FEM DOFs | Adjoint wall time (s) | FD wall time (s) | Aggregate FD CPU time (s) |
|---|---:|---:|---:|---:|---:|
| pilot-128 | 128 | 13359 | 1.990 | 120.379 | 119.795 |
| scale-32 | 32 | 52347 | 10.672 | 40.339 | 202.702 |
| scale-128 | 128 | 52347 | 12.472 | 112.083 | 755.298 |
| scale-512 | 512 | 52347 | 12.741 | 395.261 | 2999.272 |

The serial pilot speedup was 60.48× on the same node. The three scaling cases used eight FD workers on one fixed analysis mesh. Each perturbation supplied both responses.

## KKT audit

The archived β = 8 design achieved normalized power 0.992919627 and scaled side-flow constraint -9.97890521. The side requirement was inactive, with 0.997891 percentage points of slack and multiplier 0.
Projected stationarity was 1.369e-08. The primal, dual and complementarity residuals were 0, -0 and 0. The design had 49 lower-bound, 73 upper-bound and 6 interior variables.
The multiplier gives zero local marginal pumping cost for a sufficiently small increase in the side-flow requirement. The KKT result applies to the immersed β = 8 problem. MMA continuation was unnecessary.
The history replayed 97 archived unique MMA evaluations.

## Mesh and hydraulic comparison

| Geometry | Target spacing (mm) | Triangles | FEM DOFs | Power (µW/m) | Resistance (Pa·s/m²) | Side flow (%) |
|---|---:|---:|---:|---:|---:|---:|
| immersed-original | 0.25 | 2848 | 13359 | 8.893614 | 3530.514 | 7.414738 |
| immersed-original | 0.125 | 11392 | 52347 | 8.892091 | 3529.910 | 7.431694 |
| immersed-original | 0.0625 | 45568 | 207219 | 8.892020 | 3529.881 | 7.437211 |
| immersed-raw | 0.25 | 2848 | 13359 | 8.830644 | 3505.517 | 8.412628 |
| immersed-raw | 0.125 | 11392 | 52347 | 8.831105 | 3505.700 | 8.424696 |
| immersed-raw | 0.0625 | 45568 | 207219 | 8.831774 | 3505.965 | 8.429885 |
| fitted-original | 0.25 | 4132 | 19137 | 9.282894 | 3685.047 | 4.038784 |
| fitted-original | 0.125 | 16466 | 75180 | 9.277632 | 3682.959 | 4.027304 |
| fitted-original | 0.0625 | 65970 | 299028 | 9.275934 | 3682.284 | 4.024206 |
| fitted-raw | 0.25 | 4338 | 20154 | 8.834575 | 3507.077 | 8.443024 |
| fitted-raw | 0.125 | 16986 | 77570 | 8.828640 | 3504.721 | 8.459634 |
| fitted-raw | 0.0625 | 67860 | 307548 | 8.827238 | 3504.165 | 8.464686 |
| fitted-smoothed | 0.25 | 4236 | 19625 | 8.867440 | 3520.124 | 8.196856 |
| fitted-smoothed | 0.125 | 16961 | 77435 | 8.859554 | 3516.993 | 8.221901 |
| fitted-smoothed | 0.0625 | 67865 | 307533 | 8.858017 | 3516.383 | 8.229312 |

| Comparison at finest paired level | Power saving (µW/m) | Saving (%) | Side gain (pp) | Power uncertainty screen (µW/m) | Raw energy residual sum (µW/m) | Mesh screen |
|---|---:|---:|---:|---:|---:|---|
| immersed-raw | 0.060246 | 0.6775 | 0.992674 | 0.001984 | 0.002994 | Pass |
| fitted-raw | 0.448696 | 4.8372 | 4.440481 | 0.011197 | 0.001788 | Pass |
| fitted-smoothed | 0.417917 | 4.5054 | 4.205106 | 0.013148 | 0.001788 | Fail |

The power uncertainty screen uses successive paired and individual changes. It is a descriptive discretization estimate. The finite-resistance immersed model has a minimum nominal penetration length of 0.10010 mm and the 0.0625 mm target provides about 1.60 nominal intervals across that length.

The smoothed wall's coarse-to-medium side-flow change was 0.02505 percentage points against a 0.020 target. Its result remains provisional under the two-comparison criterion. The primary raw fitted geometry met the three-level screen, so the secondary smoothed geometry was left at the planned three levels.

The fitted areas are 83.000, 85.507 and 85.458 mm² for the original, raw and smoothed geometries. The corresponding maximum wall turns over 0.25 mm are 90.000°, 54.678° and 55.558°. Area and outlet share both change, so the study cannot isolate the effect of turning angle. No compatible matched-area and matched-flow control was available.

## Reproducing the evidence

The complete scripts, frozen inputs, machine-readable results, PDF and PNG figures, Slurm logs and recorded job identifiers are included. See the method notes next to the verification scripts. The staging launcher reproduces the main Slurm studies when run in a fresh copy of the staged input and source tree. Run `bash rebuild_figures.sh` to regenerate the plots and tables with the pilot Python dependencies available.

All 28 common pilot input files matched the supplied Fixed_Pipe_Pilot.zip byte for byte. The archive fingerprints and the diffs for the three modified source files are included. The changes expose independent mesh and density-grid spacing and allow response-only evaluation for finite differences. The original defaults reproduced the archived mesh and responses. The topology and physical smoothed contour were frozen for these studies.
