# Coolant Manifold — Proposed Sensitivities and Control Parameters

## 1. Research Objective

Investigate how Brinkman density and other design parameters influence the hydraulic performance and reconstructed wall geometry of a two-dimensional coolant manifold.

The initial study will focus on three sensitivities:

1. Side-port flow-rate sensitivity.
2. Pumping-power sensitivity.
3. Reconstructed wall-position sensitivity.

The eventual goal is to use these sensitivities in a gradient-based topology optimization framework coupled to a validated flow solver.

---

## 2. Proposed Sensitivities

### 2.1 Side-Port Flow-Rate Sensitivity

**Quantity of interest:**

\[
J_1 = Q_{\mathrm{port}}
\]

**Sensitivity:**

\[
G_1 = \frac{dQ_{\mathrm{port}}}{d\boldsymbol{\gamma}}
\]

This measures how changes in the Brinkman density field affect flow through the perpendicular side port.

**Research questions:**

- Which regions most strongly influence the flow split?
- Where would increasing fluid density improve side-port flow?
- How do the sensitivity patterns change with operating conditions?

**Proposed method:** Discrete adjoint sensitivity analysis, verified against finite differences.

### 2.2 Pumping-Power Sensitivity

**Quantity of interest:**

\[
J_2 = P_{\mathrm{pump}}
\]

**Sensitivity:**

\[
G_2 = \frac{dP_{\mathrm{pump}}}{d\boldsymbol{\gamma}}
\]

This measures how design changes affect the hydraulic power required to maintain the prescribed inlet flow.

For a single inlet and multiple outlets, a useful pressure-power measure is:

\[
P_{\mathrm{pump}} =
p_{\mathrm{in}}Q_{\mathrm{in}}
-
\sum_{k=1}^{N_{\mathrm{out}}}p_kQ_k
\]

Here the pressures must be defined consistently as appropriate boundary or flow-weighted values.

**Research questions:**

- Which regions contribute most to pumping-power requirements?
- Does increasing side-port flow also increase hydraulic resistance?
- Can the manifold achieve a desired flow split without excessive pumping power?

**Proposed method:** Discrete adjoint sensitivity analysis, verified against finite differences.

### 2.3 Reconstructed Wall-Position Sensitivity

**Quantity of interest:**

\[
J_3 = y_w
\]

**Local geometric sensitivity:**

\[
G_3 = \frac{\partial y_w}{\partial\boldsymbol{\gamma}}
\]

This measures how changes in Brinkman density affect the location of the reconstructed solid–fluid boundary.

The proposed wall-extraction method includes:

1. Identify intermediate-density cells associated with a solid–fluid transition.
2. Calculate an average density threshold for that transition.
3. Locate the threshold crossing using linear or quadratic interpolation.
4. Reconstruct the wall and calculate the wall-distance field.

For linear interpolation:

\[
y_w =
y_1+
\frac{\bar{\gamma}-\gamma_1}
{\gamma_2-\gamma_1}(y_2-y_1)
\]

**Research questions:**

- How sensitive is wall location to changes in density?
- Does quadratic interpolation produce a more stable reconstructed wall than linear interpolation?
- Does column-average thresholding create discontinuities in the geometric sensitivity?
- How do reconstructed wall sensitivities compare with those obtained using a fixed density threshold?

**Proposed method:** Direct differentiation of the wall-reconstruction formula, verified using finite differences.

Wall position is a geometric response rather than a flow objective. A scalar wall-error measure could later be introduced if an adjoint-based flow or geometry optimization requires it.

**Important:** The wall-extraction method must handle multiple interfaces within a column and nonvertical interfaces near the side-port junction.

---

## 3. Potential Controlling Parameters

Not every adjustable parameter should be optimized. Some define the physical design, while others specify operating conditions or numerical approximations.

| Parameter | Symbol | Function | Proposed role |
|---|---|---|---|
| Brinkman density field | \(\gamma(x,y)\) | Controls fluid/solid distribution | Primary design variable |
| Side-port position | \(x_s\) | Changes branch location | Future design variable |
| Side-port width | \(w_s\) | Changes branch opening | Future design variable |
| Channel height | \(H\) | Changes available flow area | Optional geometric variable |
| Channel length | \(L\) | Changes manifold dimensions | Optional geometric variable |
| Inlet velocity | \(U_{\mathrm{in}}\) | Controls incoming flow | Operating parameter |
| Main outlet pressure | \(p_m\) | Controls downstream pressure | Operating parameter |
| Side-port pressure | \(p_s\) | Controls side-port pressure | Operating parameter |
| Maximum Brinkman resistance | \(\alpha_{\max}\) | Controls solid-region penalization | Numerical/model parameter |
| Penalization exponent | \(q\) | Controls intermediate-density resistance | Numerical/model parameter |
| Density filter radius | \(r_f\) | Controls design smoothing and feature scale | Numerical/design regularization |
| Projection sharpness | \(\beta\) | Controls density projection | Numerical continuation parameter |
| Wall-extraction threshold | \(\bar{\gamma}\) or \(\gamma_c\) | Determines reconstructed wall position | Method-comparison parameter |
| Turbulence-model coefficients | Model-dependent | Control modeled turbulence behavior | Model-validation parameters |

### 3.1 Primary Design Controls

Initially, the design vector should consist only of the Brinkman density field:

\[
\mathbf c =
[\gamma_1,\gamma_2,\ldots,\gamma_N]^T
\]

Each component controls the density of a design cell.

### 3.2 Future Geometric Controls

A later study could expand the design vector:

\[
\mathbf c =
[\gamma_1,\ldots,\gamma_N,x_s,w_s]^T
\]

This would allow the optimizer to modify the internal topology and the side-port geometry.

### 3.3 Operating Conditions

Inlet velocity and outlet pressures should initially remain fixed.

Separate parameter studies can investigate how operating conditions influence the sensitivity fields and optimal geometry.

### 3.4 Numerical and Modeling Controls

Parameters such as \(\alpha_{\max}\), \(q\), \(r_f\), and \(\beta\) should be investigated for robustness and convergence.

They should not initially be optimized as physical design variables because doing so could improve a numerical objective without producing a better physical manifold.

---

## 4. Adjoint Sensitivity Framework

Represent the discretized governing equations as:

\[
\mathbf R(\mathbf z,\mathbf c)=0
\]

where \(\mathbf z\) contains the flow solution variables and \(\mathbf c\) contains the selected controls.

For a scalar response \(J\), solve the adjoint equation:

\[
\left(\frac{\partial\mathbf R}{\partial\mathbf z}\right)^T
\boldsymbol\lambda
=
\left(\frac{\partial J}{\partial\mathbf z}\right)^T
\]

The total response sensitivity is:

\[
\frac{dJ}{d\mathbf c}
=
\frac{\partial J}{\partial\mathbf c}
-
\boldsymbol\lambda^T
\frac{\partial\mathbf R}{\partial\mathbf c}
\]

Separate adjoint solutions will generally be required for side-port flow and pumping power.

Wall-location sensitivities can be calculated directly from the reconstruction algorithm. Their effect on flow responses must also be included when the reconstructed wall participates in the governing equations.

---

## 5. Proposed Implementation Sequence

### Phase 1 — Fixed-Geometry Flow Verification

- Complete the fixed-geometry 2D manifold solver.
- Verify mass conservation, pressure conditions, and outlet flow calculations.
- Verify pumping-power calculations.
- Establish mesh-convergence behavior.

### Phase 2 — Brinkman Sensitivities

- Introduce a prescribed Brinkman density field.
- Calculate side-port flow and pumping-power sensitivities.
- Compare tangent, adjoint, and finite-difference results.
- Perform directional Taylor-remainder tests.

### Phase 3 — Wall-Reconstruction Sensitivities

- Implement transition-specific column-average thresholding.
- Implement linear and quadratic interpolation.
- Compare against fixed-threshold extraction.
- Calculate derivatives of reconstructed wall position.
- Check continuity, root selection, and gradient accuracy.

### Phase 4 — Coupled Wall and Flow Sensitivities

- Incorporate reconstructed walls and their boundary conditions into the flow model.
- Account for wall-position and wall-distance dependence in sensitivity calculations.
- Validate against explicit-wall simulations.
- Assess mesh and interpolation sensitivity.

### Phase 5 — Optimization

- Introduce a constrained optimization algorithm.
- Optimize the Brinkman density field.
- Investigate flow split versus pumping-power tradeoffs.
- Consider side-port position and width as additional controls.

---

## 6. Initial Research Recommendation

Prioritize the following three response gradients:

\[
\boxed{
\frac{dQ_{\mathrm{port}}}{d\boldsymbol\gamma},
\qquad
\frac{dP_{\mathrm{pump}}}{d\boldsymbol\gamma},
\qquad
\frac{\partial y_w}{\partial\boldsymbol\gamma}
}
\]

Keep geometric dimensions and operating conditions fixed during initial verification.

Only introduce additional design controls after the flow solver, adjoint sensitivities, and wall reconstruction have been independently validated.

**Status:** Proposed research plan. This document does not claim that the complete RANS–Brinkman optimization or wall-sensitivity framework has already been implemented.
