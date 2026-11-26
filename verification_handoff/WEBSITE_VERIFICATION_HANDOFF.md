# Website verification handoff
Use raw beta8 Re50 NS candidate; planar CFD, cylindrical display only. Power µW/m depth, resistance P/Qin², side share Qside/Qin. Finest fitted original→raw: 9.275934→8.827238 µW/m, 3682.284→3504.165 Pa s/m², 4.024206→8.464686% side; 4.8372% saving and +4.440481 pp. Area83→85.507, turn90→54.678°. Immersed optimizer reduction .708037%, fine .6775%, not fitted saving. Smoothed 8.858017,8.229312%, area85.458 turn55.558°, failed .02505 pp vs .020 mesh gate. Provisional. 
Include audit ladder, actual/exact/error/gate/status, historical source caveat; all 15 mesh rows and mass/energy/leakage metrics in JSON. KKT beta0/2/8 17/18/62 unique; 97 total; f=.992919627, g=-9.97890521, lambda0, slack .997891pp, stationarity1.369e-8, bounds49/73/6. Final only one PDE inequality, no area/turn cap. Adjoint 15 FD passes, original Taylor9/15, follow-up6/6 orders2.0010–2.0060, full coordinate errors2.376e-12/1.099e-12, serial60.48×; separate parallel/CPU work. Stokes/NS all six cross-evaluations and fitted three-level table in [nonlinear report](../evidence/nonlinear/nonlinear_component.md). Stokes under NS fails side floor; finest fitted Stokes9.105665µW/m and1.948290% versus NS8.827238 and8.464686%.
## Existing figures — exact paths
| Purpose | Path |
|---|---|
|Taylor initial|evidence/verification/adjoint-evidence/adjoint_taylor_remainders.png|
|Taylor follow-up|evidence/verification/adjoint-coordinate-evidence/adjoint_taylor_remainders.png|
|Sensitivity cost|evidence/verification/timing-2408534/figures/fixed_pipe_sensitivity_cost.png|
|KKT history|evidence/verification/kkt-results/fixed_pipe_kkt_history.png|
|KKT spatial|evidence/verification/kkt-results/fixed_pipe_kkt_spatial.png|
|Mesh convergence|evidence/verification/mesh-evidence/mesh_convergence.png|
|Fitted pressure/speed|evidence/verification/mesh-evidence/fitted_pressure_speed.png|
|Immersed pressure/speed|evidence/verification/mesh-evidence/immersed_pressure_speed.png|
|Junction|evidence/nonlinear/figures/junction_comparison.png|
|Wall overlay|evidence/nonlinear/figures/optimized_walls_overlay.png|
|Density difference|evidence/nonlinear/figures/density_difference.png|
|Velocity|evidence/nonlinear/figures/velocity_comparison.png|
|Pressure|evidence/nonlinear/figures/pressure_comparison.png|
|Objective history|evidence/nonlinear/figures/objective_history.png|
|Side history|evidence/nonlinear/figures/side_flow_history.png|
|Runtime|evidence/nonlinear/figures/optimization_runtime.png|
Render equations from audit sections1–5: PDE/weak form/full convective derivative, anchored filter/projection/rational drag and their derivatives, port power and side quotient, transpose adjoint and projected KKT. State limitations and actual AI disclosure. No global optimum, fixed-area causal comparison, 3D CFD validation, or SST convergence. Website agent owns viewer rubric/reproducibility. Do not alter numerical archive.
