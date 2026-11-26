---
layout: default
title: "PDE-Constrained Design Optimization of a Coolant Flow Manifold"
---

# PDE-Constrained Design Optimization of a Coolant Flow Manifold

**Arizona State University — Design Optimization, Project 4**  
**Team:** [TODO: team name and authors]  
**Repository:** [TODO: public repository URL]  
**Report website:** [TODO: GitHub Pages URL]  
**Draft evidence cutoff:** November 25, 2026, latest reported width-constrained experiments not yet complete.

> **Draft status:** This report separates archived benchmark data, recently reported Sol results and pending validation. It does not assert that a physically feasible topology optimum has been obtained. Replace every TODO with actual evidence, not a predicted value. The latest width-only run is the designated follow-up; earlier experiments are retained as comparisons, not silently overwritten.

## 1. Executive summary

Coolant manifolds distribute fluid between channels while consuming mechanical energy supplied by a pump. This study investigates whether a gradient-based optimizer can reshape the fluid–solid boundary of a planar T-junction through a Brinkman density field, rather than prescribing individual fillet radii or baffles. The engineering aim is to reduce hydraulic losses, preserve a usable side-outlet flow and avoid excessive passage restrictions or sharp wall transitions.

The team developed a custom mixed finite-element solver for steady incompressible Stokes and Navier–Stokes flow, implemented a discrete adjoint and used gradient-based optimization with density filtering and projection. The nonlinear study uses an inlet Reynolds number of 50. The reference manifold has a 20 mm by 4 mm main passage and a 1 mm by 3 mm branch. The original physical fluid area is 83 mm².

The reported 1 mm wall-editing study contained 960 density entries, of which 416 were independent design variables. Eight independent derivative cases passed for each of the linear and nonlinear solvers. The nonlinear run reached an independently KKT-stationary endpoint at $\beta=4$, but its native evaluation-limit exit blocked advancement under the then-current execution policy. The reconstructed geometry remained connected and reduced the maximum measured wall turn from 90° to 66.44°. Its fitted-wall flow fraction, however, fell below the prescribed 10% minimum, and its physical pumping power was not lower than the original sharp reference. A feasible Stokes endpoint was not established. A width-only follow-up isolates the effect of prohibiting substantial contractions without requiring expensive geometric finite differences.

| Headline quantity | Original fitted sharp T | Diagnostic $\beta=4$ immersed endpoint | Same endpoint fitted replay |
|---|---:|---:|---:|
| Mechanical port power (µW/m) | 9.27593 | 8.7852 | 9.3114 |
| Side-outlet flow (%) | ≈4.027 | 10.0000 | 9.3628 |
| Physical area (mm²) | 83 | 87.3897 (reconstructed) | 87.3897 |
| Maximum wall turn at 0.25 mm scale | 90° | 66.44° (reconstructed) | 66.44° |
| Physical engineering certification | Reference, not feasible for 10% target | Not applicable to porous state | Not achieved |

The immersed endpoint is 5.29% below the sharp fitted power numerically, whereas the fitted replay is 0.38% above it. Neither is a matched-performance savings claim. The two geometries have different flow splits and fluid areas. The purpose of the subsequent experiment is to resolve whether a narrower design space can produce a more useful geometry, not to retroactively declare the earlier result successful.

### Interactive 3D comparison (required at the top of the report)

<div id="manifold-viewer" style="min-height:420px" aria-label="Initial and candidate manifold comparison">
<p><strong>TODO:</strong> Embed the live WebGL comparison here using the exported solver data. Show the original geometry and the most recent candidate, with honest labels for incomplete or failed constraints.</p>
</div>

<!-- TODO: connect the actual viewer module after it exists; do not publish a broken script link. -->

Display the 2D analysis geometry extruded by an explicitly identified visualization depth; extrusion is not a new three-dimensional flow simulation. Use a common field scale, show the flow and power metrics, and mark failures independently of color. See Section 11 for the complete implementation and submission checklist.

## 2. Engineering structure and optimization formulation

### 2.1 Motivation and scope

The manifold is a representative single-inlet/two-outlet coolant distribution component. A sharp junction may introduce strong velocity gradients and local recirculation. The channel dimensions and the relative resistance of both outlet pathways also determine flow allocation. Changing the wall shape can therefore influence both mechanical energy loss and side-outlet delivery. Smoothness alone is not a guarantee of lower power.

The current model is a planar hydraulic study, not a full conjugate heat-transfer or turbulent cooling simulation. No thermal-performance claim is made. The project uses the course's listed flow-manifold option and its N2 nonlinear physics component.

### 2.2 Geometry and fixed design regions

Coordinates below use millimeters. The reference fluid region is

$$
\Omega_0=([0,20]\times[0,4])\cup([9.5,10.5]\times[4,7]).
$$

The inlet is at $x=0$, the main outlet is at $x=20$, and the side outlet is at $y=7$ with $9.5\le x\le10.5$. The junction corners are $(9.5,4)$ and $(10.5,4)$. The design rectangle is $\Omega_d=[6,16]\times[0,6]$, with area 60 mm²; its union with the reference passage has area 101 mm². The original design rectangle contains 42 mm² fluid and 18 mm² solid. The remaining 41 mm² of fluid is fixed.

For cell size $h_d=0.25$ mm, there are $N=40\times24=960$ density entries. The indicator $m_i=1$ represents fluid, and $m_i=0$ represents solid. Intermediate values describe a penalized porous numerical material, not a physical mixture prescribed for manufacture.

| Design quantity | Value/bounds | Rationale |
|---|---|---|
| Raw density $m_i$ | Dimensionless, $0\le m_i\le1$ | Numerical fluid indicator |
| Nominal editable window | $[6.5,15.5]\times[0,5.5]$ mm | Retains fixed collars at design interfaces |
| Wall-edit distance $d_{\max}$ | 1 mm | Team-selected geometric preservation limit |
| Original protected core | $[0,20]\times[1,3]$ plus branch core | Earlier wall-displacement run |
| Width-only protected core | $[0,20]\times[0,3.75]$ plus branch core | Allows at most 0.25 mm reduction of the original main width while fixing its bottom |
| Branch core | $[9.5,10.5]\times[2,7]$ mm | Preserves original branch connectivity and opening |
| Minimum side flow / maximum side flow | 10% / 90% | Team-selected investigation target, not an external standard |
| Physical area | Measured, no cap in latest experiments | Area cap removed to focus on wall geometry |
| Turning-angle target | 45° at 0.25 mm arc-length scale | Team-selected smoothness diagnostic; hard constraint only in job 2407390 |

The distance to a wall segment $[\mathbf a,\mathbf b]$ is computed using

$$
t=\operatorname{clip}\left(\frac{(\mathbf x-\mathbf a)\cdot(\mathbf b-\mathbf a)}{\|\mathbf b-\mathbf a\|^2},0,1\right),\qquad d_j=\|\mathbf x-\mathbf a-t(\mathbf b-\mathbf a)\|.
$$

The admissible editing region is the intersection of the wall-distance band with the editable window, excluding the protected core. A conservative implementation permits a cell only when all its vertices are within the distance limit of a common original wall segment and the cell does not intersect the protected core interior. All noneditable entries are fixed at the original binary value. This controls where modifications can occur, rather than imposing a two-sided Hausdorff-distance condition or guaranteed curvature bound.

For the current width-only experiment, the protected strip must be included in the reconstructed geometry. This provides a minimum contiguous main-channel width of 3.75 mm without introducing a separate differentiable width response. The 10% lower flow target may be difficult to reach after restricting contractions. The result will decide feasibility.

### 2.3 Objective: mechanical port power

Let the outward unit normal on each port be $\mathbf n$. Define signed inlet and side-outlet flows per unit depth by

$$
Q_{in}=-\int_{\Gamma_{in}}\mathbf u\cdot\mathbf n\,ds,\qquad Q_s=\int_{\Gamma_s}\mathbf u\cdot\mathbf n\,ds.
$$

The mechanical port-power functional (positive for a net input of pressure and kinetic energy) is

$$
P=-\int_{\Gamma_{ports}}\left(p+\frac12\rho |\mathbf u|^2\right)(\mathbf u\cdot\mathbf n)\,ds,\qquad f=P/P_{ref}.
$$

Here $P$ is in W/m, not W for a specified three-dimensional device. The fixed reference used in the recent formulation is the original fitted sharp power, approximately 9.276 µW/m. TODO: confirm the exact normalization constant from each run's configuration; do not infer it from a rounded plot.

At prescribed inlet flow, the effective energy-based hydraulic resistance $R_{eff}=P/Q_{in}^2$ has the same minimizer as $P$; it is a reported diagnostic, not another objective. The actual mechanical port-power functional is not in general equal to just one static pressure drop multiplied by inlet flow because there are two exits and kinetic fluxes.

### 2.4 Current width-only mathematical problem

$$
\begin{aligned}
\min_{\mathbf m_E}\quad & f(\mathbf z(\mathbf m_E))\\
\text{subject to}\quad &R(\mathbf z,\mathbf m_E)=0,\\
&g_1=(0.10-\eta)/0.001\le0,\\
&g_2=(\eta-0.90)/0.001\le0,\\
&\eta=Q_s/Q_{in},\\
&0\le m_i\le1,\quad i\in E,\\
&m_i=m_{0,i},\quad i\notin E.
\end{aligned}
$$

The fixed-cell mask encodes the width and wall-edit constraints. Connectivity and physical flow feasibility are checked independently after threshold reconstruction. The turning angle is a diagnostic in this width-only experiment. The earlier job 2407390 adds a smoothness inequality as described in Section 6. The corner-only expansion experiment was proposed but put on hold; no results from it are included.

### 2.5 Distinguish experiment versions

| Version | Area cap | Side-flow target | Geometric restrictions | Evidence |
|---|---|---|---|---|
| Earlier density optimization | 83 mm² | Earlier strict or relaxed bands | Broader topology freedom | Historical study; not merged with current results |
| Wall editing | None | 10–90% | 1 mm wall band, original protected core | Jobs 2406318 (NS), 2406647 (Stokes) |
| Width + turn | None | 10–90% | 3.75 mm core, 45° aggregate turn constraint | Job 2407390; Phase I ongoing at last update |
| Width only | None | 10–90% | 3.75 mm core; turn measured only | TODO: run ID and results |
| Corner-only expansion | None | Proposed baseline-preserving flow | Only local additions near corners | Held; not executed in available evidence |

The separate experiments answer different questions and should not be pooled as one converged optimization.

## 3. Analysis model

### 3.1 Governing equations and constitutive law

Steady incompressible Newtonian flow is modeled on the fixed computational region $\Omega=\Omega_0\cup\Omega_d$ during immersed optimization:

$$
\rho(\mathbf u\cdot\nabla)\mathbf u-\nabla\cdot(2\mu D(\mathbf u))+\nabla p+\alpha(\bar\gamma)\mathbf u=0,\qquad \nabla\cdot\mathbf u=0,
$$

$$
D(\mathbf u)=\tfrac12(\nabla\mathbf u+\nabla\mathbf u^T),\qquad \sigma=-pI+2\mu D(\mathbf u).
$$

Removing the convective term gives the Stokes comparison model. The fitted-wall solver instead uses only the reconstructed fluid polygon and no artificial solid drag.

| Physical input | Value |
|---|---:|
| Density $\rho$ | 998.2 kg/m³ |
| Dynamic viscosity $\mu$ | 0.001002 Pa·s |
| Characteristic length $H_c$ | 0.004 m |
| Reynolds number $Re=\rho U H_c/\mu$ | 50 for NS |
| Mean inlet speed $U$ | 0.01254758565 m/s |
| Inlet volume flux per depth $UH_c$ | 5.01903426 × 10⁻⁵ m²/s |
| Pressure scale $\mu U/H_c$ | 0.0031431702 Pa |
| Power scale per depth $\mu U^2$ | approximately 1.578 × 10⁻⁷ W/m |
| Pilot $\alpha_{\max}$ | 10⁵ kg/(m³·s) |

Use the inlet profile $u_x(y)=6U(y/H_c)(1-y/H_c)$ and $u_y=0$ at the inlet. No-slip conditions apply on impermeable external walls. Both outlet planes use the natural zero symmetric traction condition $\sigma\mathbf n=0$. This is not interchangeable with prescribing zero static pressure and zero normal velocity gradient. There is no imposed downstream resistance curve or body force in the manifold study.

With $\mathbf x=H_c\mathbf x^*$, $\mathbf u=U\mathbf u^*$, and $p=(\mu U/H_c)p^*$, the dimensionless Brinkman coefficient is $\alpha^*=H_c^2\alpha/\mu$. The dimensionless inertia coefficient is $Re$. For comparisons retain the same imposed dimensional inlet profile and scales, changing only the convection term; do not set the inlet speed to zero in Stokes.

### 3.2 Finite-element formulation

The custom solver uses triangular Taylor–Hood elements: continuous quadratic velocity and continuous linear pressure, $\mathbf u_h\in[P_2]^2$, $p_h\in P_1$. Six local nodes interpolate each velocity component and three vertices interpolate pressure. The pressure sign convention must match the assembled residual and adjoint. In dimensionless variables, the weak momentum and continuity equations are

$$
\begin{aligned}
R_u(\mathbf v)={}&\int_\Omega Re(\mathbf u\cdot\nabla)\mathbf u\cdot\mathbf v\,d\Omega
+\int_\Omega 2D(\mathbf u):D(\mathbf v)\,d\Omega\\
&-\int_\Omega p\nabla\cdot\mathbf v\,d\Omega
+\int_\Omega \alpha^*\mathbf u\cdot\mathbf v\,d\Omega,\\
R_p(w)={}&\int_\Omega w\nabla\cdot\mathbf u\,d\Omega.
\end{aligned}
$$

Dirichlet data are eliminated or otherwise enforced by the same consistent reduced residual in forward and adjoint solves. With $\mathbf z$ denoting the free state degrees of freedom, the assembled Newton step solves

$$
R_{\mathbf z}\Delta\mathbf z=-R,\qquad \mathbf z^{k+1}=\mathbf z^k+\omega_k\Delta\mathbf z.
$$

The convective directional tangent is

$$
D[(\mathbf u\cdot\nabla)\mathbf u][\delta\mathbf u]
=(\delta\mathbf u\cdot\nabla)\mathbf u+(\mathbf u\cdot\nabla)\delta\mathbf u.
$$

Both terms are required. Sparse assembly and direct or existing verified sparse solves are reused. The archived fitted sharp baseline at 0.25 mm had 2,656 triangles and 12,495 DOFs. TODO: insert the actual analysis mesh counts, quadrature rules, software versions, tolerances and runtimes for every final optimization mesh from its manifest, rather than equating the design grid with the FEM mesh.

### 3.3 Anchored density filtering and projection

The Helmholtz filter uses $H=M+r_f^2 K$, $r_f=0.5$ mm. For editable and fixed indices $E,F$, preserve $\tilde\gamma_F=m_{0,F}$ and solve

$$
H_{EE}\tilde\gamma_E=M_E m_E-H_{EF}m_{0,F}.
$$

This assumes the existing cell-centered diagonal mass operator. The fixed-cell values anchor the filter; unlike the unanchored no-flux operator, global density mass need not be conserved. The projection is

$$
\bar\gamma=\pi_\beta(\tilde\gamma)=\frac{\tanh(\beta/2)+\tanh[\beta(\tilde\gamma-0.5)]}{2\tanh(\beta/2)},
$$

with exact identity for $\beta=0$. Resistance is

$$
\alpha(\bar\gamma)=\alpha_{\max}\frac{q(1-\bar\gamma)}{q+\bar\gamma}.
$$

All fixed binary values remain binary after projection. The transfer of cell resistance to element quadrature is part of the discrete model. Its exact transpose is used in the adjoint. At the pilot resistance, the penetration scale $\sqrt{\mu/\alpha_{\max}}$ is approximately 0.1001 mm. The finite resistance is not an impermeable wall: discretization and porous leakage must be examined separately.

### 3.4 Geometry extraction and physical diagnostics

Reconstruct at threshold 0.5 using the existing cell-center interpolation and marching-squares contour rule. Let $\Omega_E$ be the editable-cell union and $\Omega_r$ the threshold polygon. Preserve the original domain outside that union:

$$
\Omega_f=(\Omega_0\setminus\Omega_E)\cup(\Omega_r\cap\Omega_E).
$$

The exact edit restriction is $\Omega_f\triangle\Omega_0\subseteq\Omega_E$. Verify the protected core is fully contained, all three ports belong to the same connected component, and there are no unintended pockets or islands. Clipping can create corners; report them rather than silently smoothing them.

Physical area is $A_f=|\Omega_f|$, with $A_+=|\Omega_f\setminus\Omega_0|$ and $A_-=|\Omega_0\setminus\Omega_f|$; $A_f=83+A_+-A_-$. It is not equivalent to the sum of raw or projected densities.

Viscous and Brinkman dissipation are

$$
\Phi_v=\int_\Omega2\mu D:D\,d\Omega,\qquad \Phi_B=\int_\Omega\alpha|\mathbf u|^2\,d\Omega.
$$

For a divergence-free continuum solution, the mechanical port power equals $\Phi_v+\Phi_B$ minus viscous port traction work. In the discrete case, include the divergence defect $-\int (p+\rho|\mathbf u|^2/2)\nabla\cdot\mathbf u\,d\Omega$ and report the residual; do not presume pointwise incompressibility of Taylor–Hood velocity.

### 3.5 Assumptions

Steady, incompressible, isothermal Newtonian flow; planar geometry, unit-depth quantities; no heat transfer, turbulence closure, roughness, manufacturing tolerances or external outlet loads modeled. The inlet is fully developed. Finite Brinkman resistance is a design-discovery approximation. The fitted replay is a higher-fidelity check of the same planar model, not an experimental validation of a three-dimensional cold plate.

## 4. Verification ladder

The archived `output/project43_baseline_verification.json` embedded in `Project43_searchable.txt` reports a complete earlier solver benchmark suite. The following numbers were extracted directly from it. These tests verify specific implementations and mathematical fixtures; they do not certify every later optimization or the physical permeability approximation.

| ID | Test and expected result | Computed result/error | Archived gate |
|---|---|---|---|
| V1 | Element interpolation and partition of unity: defect 0; quadrature moments exact | Interpolation, partition and gradient-sum defects 0; moment error 1.11e−16 | Pass |
| V2 | Planar Poiseuille: nondimensional exact flux 1 | 0.9999999999999918; relative flux error 8.22e−15; velocity L2 1.18e−14 | Pass |
| V3 | Rotated Poiseuille: exact flux 1; invariant solution | 0.9999999999999974; relative flux error 2.55e−15 | Pass |
| V4 | Mass conservation: net flux defect 0 | Relative mass imbalance 8.22e−15 in original pipe | Pass |
| V5 | Manufactured solution mesh refinement | Velocity orders 3.128, 3.097; pressure orders 3.226, 3.422; gates 2.6/1.7 | Pass |
| V6 | Nonlinear manufactured response and tangent | Velocity L2 0.0017012; pressure L2 0.0010254; tangent relative error 1.076e−11 | Pass |
| V7 | Uniform Brinkman duct exact flow 0.29057005475492537 | 0.2905694293009031; relative error 2.153e−6 | Pass |
| V8 | Layered Brinkman duct at fixture coefficient 1600; exact total flow 0.17337431884801597 | 0.17337449558640106; relative error 1.019e−6; solid flow fraction 0.041108 | Pass finite-penalty benchmark, not impermeability |
| V9 | 24-variable adjoint fixture; finite difference and Taylor convergence | At step 0.01, pressure and upper-half flux directional relative errors 1.148e−7 and 2.127e−6 | Pass |
| V10 | Current wall-band NS and Stokes derivative suites | Reported 8/8 each; TODO exact errors, case IDs and step sizes | Reported pass |
| V11 | Width + turn derivative suite | Both verification jobs reported passed; TODO identify exact fixtures | Reported pass |
| V12 | Width-only mask and gradient suite | TODO | Pending |

For Poiseuille flow with pressure gradient $G=-dp/dx$, the analytical velocity and flux per depth are $u(y)=G y(H-y)/(2\mu)$ and $Q=GH^3/(12\mu)$. Derive the expected solution independently of solver code. For manufactured tests use the exact archived manufactured field and forcing; TODO: copy their mathematical expressions and boundary data from the actual benchmark implementation and have team confirm independent derivation. Do not replace with a different invented manufactured field.

![Archived manufactured-solution error convergence](project4_report_assets/manufactured_convergence.png)

| Mesh resolution parameter | Velocity L2 | Pressure L2 |
|---:|---:|---:|
| 4 | 0.0525773 | 0.0378547 |
| 8 | 0.00601501 | 0.00404450 |
| 16 | 0.000702917 | 0.000377283 |

The nonlinear benchmark residual sequence was 0.0163032, 0.00235552, 1.43605e−7, 2.68584e−16. This documents the observed convergence; do not assign a measured asymptotic Newton order without a clean order calculation.

### 4.1 Physical mesh and model verification

Earlier binary Brinkman-to-fitted parity studies failed strict physical parity gates, despite successful analytical Brinkman benchmarks. Finite permeability can alter effective hydraulic resistance and branch flow. Preserve the failed cases. The conventional reference $b=w=3$ mm passed a separate finer FEM mesh series at 11.7745% side flow and 12.4993 µW/m. Its original series had a 0.05537 percentage-point split change, beyond the 0.020 pp limit. The finer series does not erase that prior failure.

| Physical gate | Criterion | Status |
|---|---|---|
| Fitted mesh split changes | Final prescribed series, ≤0.020 percentage points per successive change | Conventional reference reported pass; diagnostic candidate only one replay |
| Fitted mesh power changes | ≤0.5% per successive change | Conventional reference reported pass; candidate pending |
| Fitted flow | 10–90%, accounting for mesh allowance | Diagnostic candidate failed at 9.3628% |
| Independent OpenFOAM | Matched geometry, inlet, physics and outlet traction | Blocked: existing outlet BC differs |
| Width-only final candidate | Actual mask, mesh and flow checks | TODO |

TODO: supply three-level tables with mesh counts, power, side share, consecutive errors, residual and mass defect. Provide both failed and accepted reference series and exact archived file names. A finer independent series and its estimator must be identified explicitly.

## 5. Sensitivity analysis

### 5.1 Discrete adjoint derivation

At a converged state, $R(\mathbf z,\mathbf m)=0$. Differentiating gives

$$
R_{\mathbf z}\frac{d\mathbf z}{d\mathbf m}=-R_{\mathbf m}.
$$

For scalar response $F$, define $R_{\mathbf z}^T\lambda_F=\nabla_{\mathbf z}F$. Then

$$
\frac{dF}{d\mathbf m}=\partial_{\mathbf m}F-R_{\mathbf m}^T\lambda_F.
$$

The fixed port locations eliminate direct design dependence of power and flow at a fixed state. The objective adjoint right-hand side follows from

$$
\delta P=-\int_{\Gamma_{ports}}[(\delta p+\rho\mathbf u\cdot\delta\mathbf u)(\mathbf u\cdot\mathbf n)+(p+\tfrac12\rho|\mathbf u|^2)(\delta\mathbf u\cdot\mathbf n)]\,ds.
$$

Divide by $P_{ref}$ for the normalized objective. Since the prescribed inlet flow is fixed,

$$
\delta\eta=Q_{in}^{-1}\int_{\Gamma_s}\delta\mathbf u\cdot\mathbf n\,ds.
$$

The constraint gradients share one flow adjoint: $\nabla g_1=-\nabla\eta/0.001$, $\nabla g_2=\nabla\eta/0.001$. Fixed-cell bounds have zero sensitivity outside the independent variable space and require no adjoint. The width constraint is imposed exactly by the protected mask, not approximated with a CFD sensitivity.

### 5.2 Density and resistance derivative

$$
\pi_\beta'=\frac{\beta\operatorname{sech}^2[\beta(\tilde\gamma-0.5)]}{2\tanh(\beta/2)},\quad
\alpha'=-\alpha_{\max}\frac{q(q+1)}{(q+\bar\gamma)^2},\quad
\frac{\partial\tilde\gamma_E}{\partial m_E}=H_{EE}^{-1}M_E.
$$

At dimensionless quadrature point $k$ the residual design action contracted with the adjoint is $s_k=-w_k\mathbf u_k\cdot\lambda_{u,k}$. With cell-to-quadrature transfer $T_E$ and diagonal projection/drag derivatives $D_\pi,D_\alpha$,

$$
\nabla_{m_E}F=\frac{H_c^2}{\mu} M_E^T H_{EE}^{-T}D_\pi D_\alpha T_E^T\mathbf s.
$$

The transpose solve uses a factorization of the fixed anchored filter. The flow and power gradients each require a transposed state solve after convergence. The constant mass and filter matrices can be reused. Do not omit the nondimensional drag scale or transpose transfer.

### 5.3 Gradient checks and Taylor tests

For a normalized feasible perturbation direction $\mathbf d$, compare $\nabla F^T\mathbf d$ to centered differences $[F(m+h\mathbf d)-F(m-h\mathbf d)]/(2h)$. At active bounds use inward directions. The first-order Taylor remainder is $|F(m+h\mathbf d)-F(m)-h\nabla F^T\mathbf d|$; expect quadratic decay in smooth regimes until cancellation.

The following actual archived table is for the earlier 24-variable benchmark measuring pressure loss and upper-half flux, **not** the current 416-variable mechanical power objective. It supports the earlier infrastructure; add the newer full-response checks separately.

| Step | Pressure adjoint | Pressure FD | Rel. error | Flux adjoint | Flux FD | Rel. error |
|---:|---:|---:|---:|---:|---:|---:|
| 0.04 | −5.15473908 | −5.15472961 | 1.837e−6 | −0.00784162409 | −0.00784135718 | 3.404e−5 |
| 0.02 | −5.15473908 | −5.15473671 | 4.592e−7 | −0.00784162409 | −0.00784155736 | 8.509e−6 |
| 0.01 | −5.15473908 | −5.15473849 | 1.148e−7 | −0.00784162409 | −0.00784160741 | 2.127e−6 |

![Archived pressure and flux adjoint Taylor tests](project4_report_assets/taylor_benchmark.png)

| Current case | Responses checked | Reported independent verdict | Numerical details needed |
|---|---|---|---|
| NS wall editing | Residual design action, port power and side flow | 8/8 passed | TODO each probe, error and Taylor slope |
| Stokes wall editing | Same response types | 8/8 passed | TODO |
| Width + turn | Fluid and geometric checks | Two jobs passed | TODO geometry finite-difference convergence and contour events |
| Width only | Fluid gradients and modified mask | Pending | TODO |

The project-specific proposed Taylor slope gate is at least 1.8. Plot the actual power and side-flow Taylor remainders for the final design map, not just the earlier benchmark. Test independent directions and active bounds.

### 5.4 Cost comparison with finite differences

For 416 editable densities, centered coordinate finite differences of both fluid responses use up to 832 forward solves (each solve can evaluate both responses); forward differences use 416 perturbation solves plus a shared reference solve. The adjoint uses one converged state and two response transpose solves, plus inexpensive algebraic/filter actions. This is a solve-count comparison, not a measured speedup. One nonlinear forward solve may itself contain several Newton linear solves.

| Timing measurement | Finite difference | Adjoint |
|---|---|---|
| Number of independent variables | 416 for original wall band | 416 |
| Full two-response centered gradient | 832 forward solves (+ reference if needed) | 1 forward + 2 transpose solves |
| Measured time on same Sol allocation | TODO seconds | TODO seconds |
| Actual speedup | TODO | TODO |

Use identical hardware, state tolerances and mesh; report number of repetitions and whether timing includes assembly, factorization and warm starts. The width-only editable count must be filled in separately. The geometry finite differences of job 2407390 do not require CFD for each perturbation but may incur substantial reconstruction time; measured timings are needed to establish the bottleneck.

## 6. Optimization method and KKT interpretation

### 6.1 Algorithm and scaling

Use the existing MMA/CCSA gradient-based constrained optimizer with a SLSQP epigraph feasibility stage. The temporary Phase I variable $t$ solves $\min t$ subject to $g_j+\delta\le t$, $t\ge0$ and the design bounds. Use the prescribed $\delta=0.05$ then zero-margin fallback. The latest prompts allow 250 unique design evaluations per attempt and 500 in Phase II; archive total evaluations too. Objective and flow residual scaling are fixed as in Section 2. The original wall-band run used the following density continuation:

| Stage | q | β | NS 2406318 |
|---:|---:|---:|---|
| 1 | 0.01 | 0 | Pass |
| 2 | 0.05 | 0 | Pass |
| 3 | 0.10 | 0 | Pass |
| 4 | 0.10 | 1 | Pass |
| 5 | 0.10 | 2 | Pass |
| 6 | 0.10 | 4 | Independent KKT pass; native limit stopped advancement |
| 7 | 0.10 | 8 | No completed result in available evidence |

The new width-only experiment starts from the archived $\beta=4$ design, enforces new bounds, and uses $q=0.10,\beta=8$ directly for a controlled comparison with the width-plus-angle experiment. It is not a repeat of the seven-stage continuation. TODO: report exact start and endpoint objective, flow, angle, mask count, iterations and elapsed time.

### 6.2 Geometric angle experiment and cost

For points a fixed arc distance $\ell=0.25$ mm before and after a wall sample, form unit secants $\mathbf t_-,\mathbf t_+$. The measured turn is $\Theta=\arccos(\mathbf t_-\cdot\mathbf t_+)$. The width-plus-angle run imposed a KS aggregate of $g_{\theta,k}=\cos45^\circ-\mathbf t_-\cdot\mathbf t_+$:

$$
G_\theta=\tau\log\sum_k\exp(g_{\theta,k}/\tau),\qquad \tau=0.001,\qquad G_\theta/0.1\le0.
$$

Use stable log-sum-exp; the unnormalized aggregate bounds the maximum sampled violation from above. This is a sampled finite-scale constraint, not a continuous curvature guarantee. The gradient in this experiment used geometry-only finite differences, including reconstruction; it has no fluid adjoint contribution. Topology changes, clipping and resampling can introduce nondifferentiability. The width-only experiment deliberately removes this nonlinear inequality and all geometric finite-difference calls, retaining the same geometric measurement after optimization.

Job 2407390 was still in the first restoration attempt after about 26 minutes, with two verification jobs passed. Starting flow 8.9040%, starting maximum turn 66.65°. No trial history had yet been flushed, and no feasibility decision was available. This is a recorded limitation of observability, not proof of an infinite loop. TODO: replace with final status, timing profile and results when supplied.

### 6.3 KKT conditions and multipliers

For $g_j\le0$, $\mathcal L=f+\sum_j\lambda_j g_j$. Check primal feasibility, nonnegative multipliers, complementarity and box-projected stationarity:

$$
r_{stat}=\|m_E-\operatorname{clip}(m_E-\nabla f-J_g^T\lambda,0,1)\|_\infty.
$$

| Independent check | Limit | NS β=4 endpoint |
|---|---:|---|
| Primal violation | 1e−6 | Reported complete audit pass; TODO numerical value |
| Dual violation | 1e−8 | TODO numerical value |
| Projected stationarity | 1e−4 | 1.27e−6, pass |
| Complementarity | 1e−5 | TODO numerical value |
| Native optimizer termination | Archive verbatim | 500 unique evaluations exhausted |
| Lower-flow multiplier $\lambda_1$ | Nonnegative | TODO actual witness |
| Upper-flow multiplier $\lambda_2$ | Nonnegative | TODO actual witness |

The lower-flow constraint is active at the reported 10% endpoint, while the upper flow constraint is slack. In the scaled problem, under local regularity, the multiplier for the lower-flow constraint gives the change in optimal normalized power per unit change in the required fraction divided by 0.001. A positive lower-bound multiplier means increasing required side flow has a local cost. The actual multiplier must be retrieved, not inferred from the active constraint. A KKT pass is a first-order result, not a proof of global optimality or physical feasibility.

### 6.4 Convergence history and engineering interpretation

TODO: insert actual objective and flow histories with continuation boundaries, feasibility markers and stage-specific normalization. Include mask and reconstructed geometry. The reported earlier candidate expanded the wall near the branch, contracted the main channel downstream and created a secondary expansion. This is consistent with the lower flow bound encouraging redistribution, though the causal contributions require matched comparisons. At fixed $Q$, a planar fully developed segment has $R_h=12\mu L/w^3$; restricting its width to at least 3.75 mm limits the ideal same-length resistance ratio relative to 4 mm to $(4/3.75)^3\approx1.214$. This local estimate is not a prediction for the entire junction.

## 7. Nonlinear component (N2): Stokes versus Navier–Stokes

The nonlinear term $\rho(\mathbf u\cdot\nabla)\mathbf u$ is the designated N2 physics. Stokes omits it, leaving a linear state equation for fixed density. The density-to-resistance map and constrained optimization remain nonlinear even when the state equation is linear. The NS tangent and adjoint are evaluated at the converged state, with both convection tangent terms retained.

| Comparison item | Stokes 2406647 | NS 2406318 |
|---|---|---|
| Same original 1 mm mask | 416 editable cells | 416 editable cells |
| Independent derivative cases | 8/8 passed | 8/8 passed |
| Flow restoration | Both attempts reached 9.96346%; below 10% | Succeeded |
| Optimized endpoint | Not obtained; MMA did not start | Diagnostic β=4 endpoint, independently stationary |
| Physical fitted replay | Not started | 9.3628%, 9.3114 µW/m at 0.125 mm |
| Matched optimum comparison | Not established | Not established |

The Stokes restoration shortfall is 0.03654 percentage points; a local restoration outcome does not prove global infeasibility. The course asks for the effect of nonlinearity on the optimized design. This remains incomplete because there is no feasible Stokes optimum. If feasible, finish the same-model optimization on both state equations. Otherwise present cross-evaluation of identical archived fields and explicitly explain the missing optimized-design comparison. Do not call an infeasible restoration design an optimum.

| Same fixed density field | Stokes power / split | NS power / split | Gradient angle / relative timing |
|---|---|---|---|
| Original raw field, same q/β | TODO | TODO | TODO |
| NS β=4 endpoint | TODO | 8.7852 µW/m / 10.0000% | TODO |
| Width-only endpoint | TODO | TODO | TODO |

Compute directional gradient agreement using the cosine $a^Tb/(\|a\|\|b\|)$ for nonzero gradients and the same independent coordinates. Evaluate both designs under nonlinear fitted physics if endpoints exist. Explain changes in flow redistribution, losses, state-solve costs and optimal geometry.

## 8. Results and physical design validation

### 8.1 Initial and candidate comparison

| Metric | Requirement | Sharp fitted reference | NS β=4 diagnostic candidate | Width-only candidate |
|---|---|---:|---:|---|
| Mechanical power (µW/m) | Minimize | 9.27593 | 9.3114 fitted; 8.7852 immersed | TODO |
| Side flow (%) | 10–90 for this study | ≈4.027, fail | 9.3628 fitted, fail | TODO |
| Main width (mm) | 3.75 for width-only | 4, pass | TODO; older mask only 2 mm core | TODO |
| Wall edit distance (mm) | ≤1 | 0 | TODO audited maximum/set inclusion | TODO |
| Connectivity | All ports connected | Yes | Yes | TODO |
| Fluid pockets / solid islands | None | None | None reported | TODO |
| Turning angle (degrees) | Diagnostic target 45 | 90 | 66.44, misses target | TODO |
| Physical fluid area (mm²) | Report only | 83 | 87.3897 | TODO |
| Physical mesh check | Prescribed finer series | TODO exact baseline series | Single replay only | TODO |
| KKT status | Complete audit | Not optimized | Passed at β=4, limit exit | TODO |

The area increase is 4.3897 mm², or 5.29%. The turning-angle decrease is 23.56° (26.18%). The side-flow fraction is higher than the original reference, but below the engineering target in the fitted replay. The immersed-to-fitted power discrepancy is 0.5262 µW/m, approximately 5.99% of the immersed value; the split discrepancy is 0.6372 percentage points. This quantifies a modeling limitation; numerical gradient accuracy does not remove it.

### 8.2 Conventional reference and fair comparison

The conventional sharp-edged comparison family expands the side width $b$ over $y=4$ to 5.5 mm and changes the main width $w$ over $x=12$ to 15.5 mm. The regions do not overlap for the proposed bounds. Its area is $A_c=83+1.5(b-1)+3.5(w-4)$ mm². The tested $b=w=3$ reference has analytic area 82.5 mm² and finest accepted fitted results of 11.7745% and 12.4993 µW/m. It is not feasible under the new 3.75 mm width requirement, so it serves as an earlier comparison, not a width-constrained matched competitor.

For a matched comparison, solve the area relation $w(b)=4+[A_f-83-1.5(b-1)]/3.5$ and match actual fitted flow by a bounded root search. Require both candidates to satisfy the relevant design constraints and report if this family cannot match the optimized area or flow. A reduced-power comparison against the conventional reference without those matches is descriptive only: the earlier candidate uses 25.50% less fitted power than this conventional reference but also delivers less side flow and has larger area. Do not call it a certified shape-only improvement.

For independent mesh screening, record at least three levels and use the maximum of the final two absolute successive power changes as a screening allowance $\epsilon$. A positive matched saving requires $P_{candidate}+\epsilon_{candidate}<P_{ref}-\epsilon_{ref}$. This allowance is an empirical indicator, not a rigorous error bound. A similar split allowance is needed for the physical flow constraint. TODO: insert actual matched comparison or explicitly state it was unavailable.

### 8.3 Final plots and conclusion

TODO: attach the actual raw/projected density fields, edit mask, before/after contours, pressure and velocity, objective/constraint histories, wall-turn locations and per-mesh validation tables. All colors, geometry labels and legends must distinguish porous and fitted results.

**Current conclusion:** The custom PDE solver and discrete adjoint have substantial benchmark evidence, and density optimization discovers nontrivial wall changes. The recent topology candidate does not establish the intended physically feasible hydraulic improvement. The width-only experiment is pending. If it succeeds, replace this conclusion with measured matched results; if it fails, retain the negative finding and discuss which constraints conflicted or which validation limits were not achieved.

## 9. Bonus work (claim only with evidence)

### Fully nonlinear Newton and implicit gradients

Potential +5: the N2 benchmark and nonlinear tangent evidence in Section 4, with converged-state adjoints in Section 5. Attach the exact nonlinear benchmark field, forcing, residual history and current nonlinear gradient tests. Eligibility is subject to instructor grading.

### Topology optimization

Potential +5: design densities represent fluid/solid; a Helmholtz filter and projection control small-scale patterns; the mask and connected-component checks restrict geometry. Earlier broad-density experiments permit material changes more generally than the current restricted-wall study. Compare against an actual fixed-topology reference, describe checkerboard controls and disclose finite permeability and threshold reconstruction errors. Do not claim this bonus solely because the final design uses a density array.

No claim for transient/history-dependent, surrogate or uncertainty bonuses. The course caps bonus credit at 15. The earlier turbulence experiments remain exploratory and are not used to claim a verified turbulent optimum.

## 10. Limitations and AI-use statement

### 10.1 Limitations

The two-dimensional model omits depth-dependent friction and three-dimensional branch effects, heat transfer and real downstream loads. Finite Brinkman drag can permit significant solid-region seepage, and sharp threshold reconstruction may introduce errors at clipping boundaries. Fitted-wall comparisons are required before engineering claims. The present turning-angle diagnostic samples a fixed physical length and does not certify a differentiable continuous surface. The current comparison uses different area and flow; no causal claim that rounding alone decreases power follows. MMA solves a nonconvex problem locally, and a native success flag, evaluation limit or independent KKT result must be distinguished from physical certification. The Stokes optimized-design comparison is incomplete. OpenFOAM uses a different outlet treatment in the available case and cannot currently serve as matched validation. Some final run reports are not available in this draft.

### 10.2 AI use and verification responsibility (team must confirm before publication)

ChatGPT assisted with mathematical formulation, discrete-adjoint derivations, geometric constraints, experimental planning and drafting the report. Codex assisted with implementing and testing the custom solver, optimization drivers, verification scripts, plotting and Sol job execution. The team set engineering constraints, evaluated the meaning of the results, reviewed mathematical derivations, and directed verification gates. TODO: team members specify exactly which equations, analytical expected values, source modules and code checks they personally derived or independently validated; do not represent AI-derived material as independently hand-derived without confirmation.

A concrete documented coding issue was found in the separate exploratory SST OpenFOAM setup: job 2357233 failed during initialization because the generated `fvSchemes` omitted the wall-distance method required by the SST model. Slurm reported completion but the CFD driver reported failure. The working copy added `wallDist { method meshWave; }`; that correction was untested in the archived integrated report. This illustrates why checking solver logs and physical output matters rather than relying on job exit codes. It is outside the primary laminar path. TODO: also identify one actual verification-caught bug in the graded FEM or adjoint code, with failing fixture, correction and passing re-test; do not invent a bug. Earlier Brinkman parity failure is a modeling failure, not necessarily a coding bug.

## 11. Website, reproduction and submission checklist

This implementation appendix operationalizes the course requirements. The report itself must be served as a public page, not just as a Markdown preview. Use the existing site layout, or a compatible static-site build with math rendering.

### Interactive viewer acceptance checklist

- [ ] Linked-camera comparison of original and latest candidate, with side-by-side panels or an explicit toggle.
- [ ] True geometry extruded for display; state the chosen visual depth and that the physics is planar.
- [ ] Field color scale with colorbar; use a sequential speed magnitude or appropriate pressure scale, and explicit pass/fail markers rather than rainbow colors.
- [ ] Show velocity or pressure on the geometry using the appropriate P2/P1 interpolation (or a documented pre-sampled export using those same shape functions).
- [ ] Animate the actual saved density/optimization history or another meaningful recorded state sequence; do not invent intermediate CFD states.
- [ ] Hover/click inspection showing element identity, pressure/velocity/density and relevant local data.
- [ ] Per-design metric panels listing objective and all constraints with explicit pass/fail/pending labels.
- [ ] Viewer reads generated `assets/results.json`, fits phone width and loads promptly.

### Reproducibility architecture (proposed, not an assertion of existing files)

```text
repository/
  index.md                         # published report
  _config.yml                     # static site settings
  _layouts/default.html           # MathJax and WebGL imports
  assets/results.json             # generated authoritative metrics and viewer data
  assets/viewer.js                # actual interactive viewer
  project4_report_assets/          # benchmark and final plots
  solver/                          # existing custom FEM and sensitivity code
  verification/                    # analytical and numerical verification fixtures
  scripts/                         # experiment, export and report generation drivers
```

TODO: record the real repository paths, source commit, environment and exact run commands. Provide a single working entry point that regenerates results.json, tables and plots; for Sol-only computations document submission and export dependencies, and test the command from a clean checkout. The final website should generate numerical displays from the same results file used by the viewer, rather than manually maintaining inconsistent copies. Draft tables here are a snapshot and must be wired to that output or generated during build. Use relative asset links under the GitHub Pages repository base path.

### Submission / completion audit

| Deliverable | Current evidence / action |
|---|---|
| Formulation and assumptions | Drafted; confirm latest run configuration |
| Six or more verification tests spanning all required types | Earlier benchmark data extracted; attach original expected fields and current case-level checks |
| Adjoint checks + Taylor plot | Earlier plot included; insert current actual power and flow checks |
| Measured adjoint vs FD runtime | TODO |
| Convergence and active multipliers | Partial; retrieve actual histories and multipliers |
| Linear vs nonlinear optimized designs | Incomplete; Stokes failed restoration |
| Final physical comparison | Incomplete; latest width-only run pending |
| Bonus evidence | Potential N2 and topology, not guaranteed |
| AI disclosure and actual graded solver bug | Team confirm and supplement |
| WebGL viewer, all eight checks | TODO implementation and smoke test |
| Single-command reproduction | TODO test |
| Public repo and live site | TODO URLs; open site privately to verify access |
| Submission | Submit both URLs through Canvas; course page lists Nov 25 as final site deadline |

## References and evidence provenance

1. ASU Design Optimization course, [Project 4: PDE-Constrained Design Optimization](https://designinformaticslab.github.io/DesignOptimization2025/project4.html), consulted November 25, 2026. Official assignment and scoring requirements.
2. Borrvall, T. & Petersson, J. (2003), *Topology optimization of fluids in Stokes flow*, International Journal for Numerical Methods in Fluids, 41, 77–107. Background reference listed on course webpage; consult original before adding paper-specific claims.
3. Gersborg-Hansen, A., Sigmund, O. & Haber, R. B. (2005), *Topology optimization of channel flow problems*, Structural and Multidisciplinary Optimization, 30, 181–192. Background reference listed on course webpage.
4. Project evidence `Project43_searchable.txt`, embedded `output/project43_baseline_verification.json` and `output/project43_baseline_run.json/report.json`: source of the archived benchmark numbers and original 0.25 mm baseline mesh metadata in this draft.
5. Project code `fem.py` and `project43_baseline.py`: custom Taylor–Hood analysis and reference geometry implementation. Check actual optimized-run code revision before attributing all details to those files.
6. `Project43_integrated_report.md`: older fitted/immersed comparisons and exploratory SST setup failure. Older parametric and turbulent findings are deliberately not combined with current topology candidate performance.
7. Reported Sol jobs 2406318, 2406647 and 2407390: later experiment status supplied in project discussion; raw latest trial histories and exact per-case derivative/KKT witnesses have not been independently reviewed in this draft. [Reported full Sol evidence archive](https://chatgpt.com/library?file_id=libfile_5f53673d8b688191a6b58798a6d772a2&account_id=personal) (private; export public-safe evidence into the repository instead of depending on this link).

**Final editorial note:** Delete this note and the draft banner only after replacing or explicitly resolving all TODOs. Keep a clear distinction between mathematical verification, immersed numerical optimization, fitted hydraulic validation, and matched-performance comparisons.
