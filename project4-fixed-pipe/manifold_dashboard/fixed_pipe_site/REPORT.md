<div class="eyebrow">ASU Design Optimization · Project 4</div>

# Fixed-pipe coolant manifold optimization

Reshaping the side-outlet transition reduced mechanical pumping power by **<span id="headline-power">…</span>** in the preliminary fitted-wall comparison. The side-outlet share rose from **<span id="headline-start">…</span> to <span id="headline-end">…</span>**. Both measurements use the same target mesh spacing. The raw reconstruction is the primary design candidate.

<div class="summary">
  <article><strong id="headline-area">…</strong><p>Fluid area of the raw reconstruction</p></article>
  <article><strong id="headline-turn">…</strong><p>Maximum wall turn of the raw reconstruction</p></article>
</div>

## Interactive 3D comparison

<section id="viewer" aria-label="Initial and optimized manifold comparison">
  <fieldset>
    <legend>Comparison controls</legend>
    <label>Geometry <select id="mode"><option value="planar">Numerical extrusion</option><option value="cylinder">Cylindrical interpretation</option></select></label>
    <label>Display <select id="field"><option>Geometry</option><option>Velocity magnitude</option><option>Pressure</option></select></label>
    <label><input type="checkbox" id="mesh"> FEM mesh</label>
    <label>Checkpoint <input id="checkpoint" type="range" min="0" max="0" step="1"></label>
    <button type="button" id="play">Play checkpoints</button>
    <button type="button" id="reset">Reset cameras</button>
    <button type="button" data-zoom="1.2" aria-label="Zoom in on both views">+</button>
    <button type="button" data-zoom="0.83333" aria-label="Zoom out on both views">−</button>
  </fieldset>
  <p id="viewer-status" class="muted" role="status">Loading geometry…</p>
  <div id="legend" class="legend" hidden><div class="legend-labels"></div><div class="legend-bar"></div></div>
  <div class="pair">
    <article class="view-panel"><h3>Original</h3><canvas tabindex="0" aria-label="Original manifold. Arrow keys orbit. Shift and arrow keys pan."></canvas><p class="view-metrics"></p><p class="inspect"></p></article>
    <article class="view-panel"><h3>Optimized</h3><canvas tabindex="0" aria-label="Optimized manifold. Arrow keys orbit. Shift and arrow keys pan."></canvas><p class="view-metrics"></p><p class="inspect"></p></article>
  </div>
</section>

Geometry mode extrudes the exact planar reconstruction by a display depth of 0.8 mm. The cylindrical interpretation wraps the measured branch half-width about the outlet axis and joins it to a 4 mm diameter main pipe. The cylindrical view represents a possible engineering geometry. Flow results come from planar simulations. Field mode shows the archived immersed β = 8 analysis plane, including the penalized part of the domain. The two views share a field scale and camera. Playback steps through actual saved continuation endpoints.

## Engineering problem and formulation

The main channel and the side branch deliver coolant to different locations. For a prescribed inlet flow, decreasing pump power while preserving side delivery can reduce operating energy. Opening the branch transition changes both local losses and the flow allocation, so the side-flow requirement belongs in the optimization.

| Region | Bounds in mm | Role |
|---|---|---|
| Main pipe | [0, 20] × [0, 4] | Fixed fluid |
| Outlet collar | [9.5, 10.5] × [6, 7] | Fixed fluid |
| Transition | [8, 12] × [4, 6] | Designable |

The transition is divided into 128 independent cells on a 0.25 mm grid. A dimensionless raw fluid indicator is bounded between zero (solid) and one (fluid). Fixed parts retain their original indicators. Let Qᵢ and Qₛ be the positive inlet and outward side flux per unit depth. The reference is the unchanged **immersed** T-junction at the same continuation stage. Its side fraction sets the allowable minimum.

\[
\min_{0\le d\le1} f=\frac{P(z(d))}{P_0},\qquad
R(z,d)=0,\qquad g=\frac{\eta_0-Q_s/Q_i}{0.001}\le0.
\]

\[
P=-\sum_{\mathrm{ports}}\int_\Gamma
\left(p+\frac{\rho}{2}|\mathbf u|^2\right)\mathbf u\cdot\mathbf n\,ds,
\quad R_{\mathrm{eff}}=\frac{P}{Q_i^2}.
\]

The normal is outward and inlet flux has the opposite orientation. Power has units W/m of out-of-plane depth. Effective resistance, area and turning angle are reported diagnostics. The pilot has one PDE-based inequality. The course guidance for team-proposed projects specifies two PDE-based constraints. The posted flow-manifold option and the pilot's implemented scope are linked below for review. [Official Project 4 assignment](https://designinformaticslab.github.io/DesignOptimization2025/project4.html)

## Governing physics and FEM

The planar model assumes steady, incompressible Newtonian flow with constant viscosity (0.001002 Pa·s) and density (998.2 kg/m³). The inlet Reynolds number based on mean speed and the 4 mm channel height is 50. The inlet has prescribed parabolic velocity. Walls have no slip. Both outlets use zero natural traction for the symmetric Newtonian stress.

\[
\rho(\mathbf u\cdot\nabla)\mathbf u-\mu\nabla^2\mathbf u+\nabla p+\alpha(\gamma)\mathbf u=0,\qquad \nabla\cdot\mathbf u=0.
\]

\[
\int_\Omega \left[2\mu\varepsilon(\mathbf u):\varepsilon(\mathbf v)
+\rho((\mathbf u\cdot\nabla)\mathbf u)\cdot\mathbf v
-p\nabla\cdot\mathbf v-q_t\nabla\cdot\mathbf u+\alpha\mathbf u\cdot\mathbf v\right]d\Omega=0.
\]

Quadratic velocity and linear pressure form a Taylor–Hood pair on triangular elements. The immersed pilot uses a structured fitted envelope with spacing at most 0.25 mm. Its design-dependent resistance is transferred cellwise. Newton solves the mixed residual with both derivatives of convection, the transported-velocity and transporting-velocity terms. For browser rendering, each scan of the β = 8 field evaluates quadratic velocity and linear pressure from the solver's local degrees of freedom. The side-flow constraint integrates the velocity trace at the port.

## Numerical verification ladder

The ladder below comes from the earlier team FEM report draft. It documents verification of common operators. The fixed-pipe audit and mesh studies have separate slots below.

<div class="table-scroll" markdown="1">

| Check | Reference or gate | Measured result |
|---|---|---|
| Element basis and integration | Partition of unity, exact polynomial moments | Defect 0, quadrature moment error 1.11 × 10⁻¹⁶ |
| Exact planar Poiseuille solution | Dimensionless flux 1 | 0.9999999999999918, relative error 8.22 × 10⁻¹⁵ |
| Rotation invariance | Rotated Poiseuille flux 1 | 0.9999999999999974, relative error 2.55 × 10⁻¹⁵ |
| Global conservation | Net flux imbalance 0 | Relative imbalance 8.22 × 10⁻¹⁵ |
| Manufactured mesh convergence | Velocity order ≥ 2.6 and pressure order ≥ 1.7 | Last velocity orders 3.128 and 3.097, pressure orders 3.226 and 3.422 |
| Nonlinear residual and tangent | Manufactured nonlinear field and FD tangent | Velocity L² error 0.0017012, pressure L² error 0.0010254, tangent relative error 1.076 × 10⁻¹¹ |
| Uniform Brinkman duct | Exact nondimensional flux 0.29057005 | Computed 0.29056943, relative error 2.153 × 10⁻⁶ |
| Layered Brinkman duct | Exact nondimensional flux 0.17337432 | Computed 0.17337450, relative error 1.019 × 10⁻⁶ |

</div>

<div class="evidence" data-evidence="immersedMesh"><strong>Immersed mesh refinement</strong><p>Expected fields: spacing, element count, power, branch share and changes under refinement.</p></div>
<div class="evidence" data-evidence="fittedMesh"><strong>Fitted-wall mesh refinement</strong><p>Expected fields: matched fixed contours, refinement levels, power and side-flow uncertainty.</p></div>

## Discrete adjoint sensitivities

For a scalar response J at a converged state, differentiate the discrete residual and eliminate the state sensitivity with a transpose solve.

\[
R_z z_d+R_d=0,\qquad R_z^T\lambda_J=J_z^T,\qquad
\frac{dJ}{dd}=J_d-\lambda_J^T R_d.
\]

For the power response f, the fixed divisor gives f_z = P_z/P₀. For the constraint, differentiate both outlet and inlet flux before applying the scale and sign.

\[
\delta P=-\sum_{\mathrm{ports}}\int_\Gamma
\left[(\delta p+\rho\,\mathbf u\cdot\delta\mathbf u)(\mathbf u\cdot\mathbf n)
+(p+\rho|\mathbf u|^2/2)\delta\mathbf u\cdot\mathbf n\right]ds,\quad
g_z=-\frac{Q_{s,z}-\eta Q_{i,z}}{0.001 Q_i}.
\]

These responses depend on d through the state. The design derivative of the residual passes through the Helmholtz filter, projection, Brinkman law and momentum assembly. If A is the free-cell Helmholtz matrix and M is its mass matrix, the filter action is A⁻¹M and the reverse action is MᵀA⁻ᵀ. The fixed-value anchors contribute only to the forward right-hand side.

\[
\gamma=\frac{\tanh(\beta\theta)+\tanh(\beta(\tilde d-\theta))}
{\tanh(\beta\theta)+\tanh(\beta(1-\theta))},\qquad
\alpha=\frac{\alpha_{\max}q(1-\gamma)}{q+\gamma},\qquad
\alpha_\gamma=-\frac{\alpha_{\max}q(1+q)}{(q+\gamma)^2}.
\]

At β = 0 the projection is the identity. The filter radius is 0.5 mm. The continuation sequence uses β = 0, 2 and 8, with q = 0.1 and maximum resistance 10⁵ kg/(m³·s). Reconstruction thresholds the projected indicator at 0.5.

At interior designs, compare each directional derivative with a central difference and measure the first-order Taylor remainder over step sizes. At active bounds use feasible directions and document which difference formula applies. One transpose solve per response gives sensitivities for all design variables. A full two-sided finite-difference gradient would require 256 perturbed state analyses for 128 variables, subject to bound feasibility. Measured timing needs a matched run.

<figure id="smoke-figure" hidden><img loading="lazy" alt="Taylor remainders for the local fixed-pipe smoke fixture"><figcaption>Local gradient smoke fixture. The independent release has its own acceptance record.</figcaption></figure>
<div class="evidence" data-evidence="finalTaylor"><strong>Final adjoint Taylor and finite-difference tests</strong><p>Expected rows: design and direction, response, steps, adjoint derivative, finite difference, error and observed order.</p></div>
<div class="evidence" data-evidence="cost"><strong>Measured adjoint versus finite-difference cost</strong><p>Expected rows: solve counts, assembly and factorization cost, wall time and speedup on matched hardware.</p></div>
<div class="evidence" data-evidence="scaling"><strong>Computational scaling</strong><p>Expected rows: degrees of freedom, design dimension, method, time and memory by mesh.</p></div>

## MMA convergence and KKT conditions

MMA builds convex moving-asymptote subproblems from the normalized power and the scaled flow constraint. The first stage starts from the original mask. Later stages warm-start from the preceding selected design. The selected checkpoint is the best feasible evaluated point identified by the pilot. [Presentation-ready optimization formulation](optimization.html)

<div id="history"></div>
<div class="table-scroll" id="checkpoints"></div>

For an inequality multiplier λ and bound multipliers μ₋ and μ₊, stationarity requires ∇f + λ∇g − μ₋ + μ₊ = 0 with primal and dual feasibility and complementary slackness. A positive λ quantifies how much the normalized objective can decrease per unit relaxation of the scaled flow inequality. The endpoint replay below checks the archived immersed endpoint and its bound stationarity. The independent audit and fitted-wall refinement have separate acceptance records.

<div class="table-scroll" id="endpoint-kkt"></div>
<div class="evidence" data-evidence="independentKkt"><strong>Independent KKT multiplier and stationarity audit</strong><p>Expected fields: binding constraints, scaled multipliers, feasibility and complementarity residuals, and engineering interpretation.</p></div>

## Nonlinear component · Stokes versus Navier–Stokes

The N2 comparison removes the convective term from the governing momentum residual to obtain Stokes flow. Use the same design space and inlet throughput and document the baseline and normalization in each solve. Evaluate each selected design under both equations to distinguish the effect of changing the geometry from the effect of changing the flow model.

<div class="evidence" data-evidence="stokes"><strong>Fixed-pipe Stokes optimization and cross-evaluation</strong><p>Expected rows: original, Stokes-selected and Navier–Stokes-selected geometries, with power, side share and constraint status under each model.</p></div>

## Fitted-wall design comparison

<figure><img src="assets/three-panel.png" loading="lazy" alt="Original sharp junction, raw reconstructed junction and smoothed junction"><figcaption>Archived comparison of the three fitted contours. The raw β = 8 reconstruction is the primary candidate.</figcaption></figure>

<div class="table-scroll" id="fitted-table"></div>

The fitted-wall study independently remeshes the three physical fluid regions at the same target spacing. The smoothing study used 10 Laplacian iterations at strength 0.5 and retained the fixed pipe and collar. Smoothness and area changed along with the outlet split. The raw shape achieved the lower power of the two optimized contours at this resolution. The fitted side-flow fraction is a separate measurement from the immersed constraint shown in the viewer.

## Bonus evidence

The nonlinear manufactured-solution test and implicit adjoint support the N2 investigation. Final-model gradient acceptance remains in the audit slot above. The topology parameterization permits fluid or solid cells within the transition and uses filtering and projection to control small features. A controlled comparison against a restricted fixed-topology optimization is still needed for the topology bonus. The fixed-pipe pilot did not investigate transient adjoints, surrogate optimization or design under uncertainty.

## Limitations and AI use

This study represents a planar laminar hydraulic model. The cylindrical visualization does not provide experimental or three-dimensional flow validation. No conclusions about heat removal, turbulence, real outlet system curves or manufacturing tolerances follow from the pilot. A finite Brinkman penalty allows residual flow in nominal solids. The fitted comparison remains preliminary until paired refinement quantifies the numerical uncertainty.

Coding agents contributed to the FEM, adjoint, MMA implementation, data exports, viewer and report. The team report draft records independent checks of analytical fluxes, manufactured solutions and quadrature. The exploratory SST setup exposed a missing wall-distance method in generated OpenFOAM input when the initialization log was checked. The team still needs to confirm the attribution and document a verified defect and regression check from the graded FEM or adjoint path.

## Reproduction and references

- [Build, numerical input manifest, Pages deployment and submission checklist](handoff.html)
- [Machine-readable numerical results](assets/results.json)
- [Official course assignment](https://designinformaticslab.github.io/DesignOptimization2025/project4.html)
- [Svanberg, The method of moving asymptotes](https://doi.org/10.1002/nme.1620240207)
- [Borrvall and Petersson, Topology optimization of fluids in Stokes flow](https://doi.org/10.1002/fld.426)
