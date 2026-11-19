# Flow Manifold Design Optimization

**Team:** AI Generated  
**Presenter:** Jack Foster  
**Project Website:** https://evl-jack.github.io/Flow.github.io/

---

## Slide 1 — Flow Manifold Design Optimization

**Team:** AI Generated  
**Presenter:** Jack Foster

**Project focus:** PDE-based fluid manifold design optimization.

---

## Slide 2 — Background

- Flow manifolds distribute fluid from a common inlet to multiple outlets.
- Internal geometry influences flow distribution, pressure loss, and pumping power.
- Poor flow distribution can reduce cooling-system performance.
- Numerical optimization offers a systematic way to improve manifold geometry.

**Visual:** Retain the existing manifold images.

---

## Slide 3 — Quick Overview

- Navier–Stokes equations with Brinkman penalization
- Fluid dynamics and computational modeling
- 2D and 3D geometry and meshing
- Numerical verification and validation
- Sensitivities, gradients, and implicit differentiation
- Gradient-based design optimization

**Workflow:**

Geometry → Mesh → Flow Solution → Sensitivities → Optimization → Validation

---

## Slide 4 — Formulation: Governing Equations

### Navier–Stokes and Brinkman Penalization

- Conservation of mass and momentum
- Brinkman resistance represents solid and fluid regions
- Design density controls local flow resistance

### Reynolds-Averaged Navier–Stokes (RANS)

- Models time-averaged turbulent flow
- Introduces Reynolds stresses through turbulence closure
- Proposed extension beyond the initial laminar model

### Shear Stress Transport (SST)

- Proposed \(k\)-\(\omega\) SST turbulence closure
- Combines near-wall \(k\)-\(\omega\) behavior with outer-flow \(k\)-\(\epsilon\) characteristics
- Requires consistent wall treatment and wall-distance calculation

---

## Slide 5 — Navier–Stokes with Brinkman

### Conservation of Mass

\[
\nabla\cdot\mathbf u=0
\]

### Conservation of Momentum

For steady, incompressible flow with constant viscosity:

\[
\rho(\mathbf u\cdot\nabla)\mathbf u
=
-\nabla p+\mu\nabla^2\mathbf u-\alpha(\gamma)\mathbf u
\]

### Brinkman Penalization

\[
\alpha(\gamma)=\alpha_{\max}(1-\gamma)^q
\]

- \(\gamma=0\): Solid-like region, high resistance
- \(\gamma=1\): Fluid region, zero Brinkman resistance
- \(0<\gamma<1\): Intermediate-density region
- \(\alpha_{\max}\): Maximum penalization strength
- \(q\): Penalization exponent

**Key concept:** The optimizer changes \(\gamma\), modifying resistance and redirecting flow.

---

## Slide 6 — Reynolds-Averaged Navier–Stokes

### RANS Momentum Equation

\[
\rho(\mathbf U\cdot\nabla)\mathbf U
=
-\nabla p+
\nabla\cdot(2\mu_{\mathrm{eff}}\mathbf S)
-\alpha(\gamma)\mathbf U
\]

Where:

\[
\mu_{\mathrm{eff}}=\mu+\mu_t
\]

\[
\mathbf S=\frac12
\left(\nabla\mathbf U+\nabla\mathbf U^T\right)
\]

### Proposed SST Turbulence Model

- Solve additional transport equations for \(k\) and \(\omega\).
- Calculate turbulent viscosity \(\mu_t\).
- Determine wall distance from the reconstructed geometry.
- Apply appropriate near-wall turbulence boundary conditions.

### Research Challenge

Brinkman resistance suppresses velocity but does not independently enforce a physical no-slip wall.

**Proposed solution:** Reconstruct the solid–fluid boundary and impose appropriate wall conditions.

**Status:** Planned RANS extension; not yet a verified production implementation.

---

## Slide 7 — Formulation: Meshing, Geometry, and Bounds

### Geometry

- Two-dimensional manifold with one inlet and multiple outlets
- Explicit polygon geometry with tagged boundaries
- Extend 2D designs into three-dimensional configurations

### Meshing

- Generate a computational mesh from the supplied geometry.
- Refine near walls, junctions, and outlet regions.
- Perform mesh-convergence studies.

### Design Bounds

\[
0\leq\gamma(x,y)\leq1
\]

- Preserve required inlet and outlet openings.
- Constrain available fluid volume or solid material.
- Apply filtering to control small geometric features.

### Proposed Wall Reconstruction

- Identify solid–fluid transition regions.
- Calculate a local intermediate-density threshold.
- Locate wall crossings using linear or quadratic interpolation.
- Compare reconstructed walls against explicit geometry.

---

## Slide 8 — Evaluation

### Mass Balance Error

\[
E_m=
\frac{|Q_{\mathrm{in}}-\sum_i Q_i|}
{Q_{\mathrm{in}}}\times100\%
\]

- Measures conservation of mass.

### Outlet Flow Uniformity

\[
Q_{\mathrm{target}}=\frac{Q_{\mathrm{in}}}{N}
\]

- Compare individual outlet flows with the equal-flow target.
- Evaluate deviations from the desired distribution.

### Pressure Drop

\[
\Delta p=p_{\mathrm{in}}-p_{\mathrm{out}}
\]

- Measures hydraulic pressure loss under the chosen pressure definition.

### Hydraulic Pumping Power

\[
P_{\mathrm{hyd}}
=
p_{\mathrm{in}}Q_{\mathrm{in}}
-\sum_i p_iQ_i
\]

- Measures pressure-power requirements.
- Enables comparison between candidate geometries.

---

## Slide 9 — Analysis & Assumptions

### Manifold Analysis Workflow

1. **Define the case:** Geometry, fluid properties, operating conditions.
2. **Build the mesh:** Identify inlet, outlets, and walls.
3. **Set up equations:** Governing physics and boundary conditions.
4. **Solve for flow:** Compute velocity and pressure fields.
5. **Calculate performance:** Flow split, pressure loss, hydraulic power.
6. **Check results:** Convergence and mass conservation.
7. **Save results:** Solution fields, metrics, verification status.

### Initial Modeling Assumptions

- Steady, incompressible flow
- Constant fluid properties
- Fixed external geometry
- Stationary no-slip walls
- Laminar baseline before turbulence-model extensions

**Visual:** Retain the existing workflow graphic.

---

## Slide 10 — Models: 2D Manifold, 1 Inlet, 3 Outlets

### Model Description

- Two-dimensional, constant-depth manifold
- One prescribed velocity inlet
- Three outlet boundaries in the presentation configuration
- Pressure or linear-resistance outlet conditions
- Explicit wall boundaries

### Governing Model

- Initial laminar, depth-averaged formulation
- Brinkman penalization planned for density-based design
- Additional plate-drag treatment for the depth-averaged model

### Current Implementation

- Polygon geometry and boundary tags defined
- Fluid-property and boundary-condition inputs implemented
- Production manifold-flow solver remains under development

**Visual placeholder:** 2D manifold geometry with labeled inlet, outlets, and walls.

---

## Slide 11 — Verification

### Numerical Verification Strategy

- Compare against analytical benchmark solutions.
- Check mesh-convergence rates.
- Verify mass, momentum, and power balances.
- Test rotational consistency.
- Verify nonlinear solver convergence.
- Compare tangent, adjoint, and finite-difference sensitivities.

### Reported Benchmark Status

- Six-stage verification program completed in follow-up testing.
- All six stages reported passing.
- Verification applies to benchmark cases, not the complete manifold solver.

**Visual placeholder:** Verification results table or mesh-convergence plot.

---

## Slide 12 — Cont.

### Verification Highlights

- Analytical strip-flow benchmarks agree closely with numerical predictions.
- Mesh refinement demonstrates expected convergence behavior.
- Conservation and rotation tests pass.
- Newton iteration demonstrates local quadratic convergence.
- Tangent and adjoint derivatives agree within numerical precision.

### Remaining Verification

- Full manifold geometry
- Multiple pressure outlets
- Brinkman density distributions
- Reconstructed-wall boundary conditions
- RANS turbulence coupling

**Key distinction:** Verified numerical components do not yet constitute a verified manifold solver.

---

## Slide 13 — Continuation / Visual

*No established title or text was present in the uploaded slide.*

**Suggested visual content:**

- Numerical verification workflow
- Benchmark geometry and solution fields
- Representative mesh-refinement comparison

**Placeholder:** Insert a figure from the numerical verification results.

---

## Slide 14 — Models: 3D Manifold

### Proposed Three-Dimensional Model

- Extend the manifold geometry into three dimensions.
- Construct explicit solid–fluid boundaries.
- Mesh internal channels and outlet junctions.
- Apply three-dimensional inlet, outlet, and wall conditions.

### Governing Physics

- Three-dimensional incompressible Navier–Stokes equations
- Proposed RANS turbulence modeling
- No-slip solid boundaries

### Purpose

- Assess three-dimensional flow effects.
- Compare 2D predictions with explicit 3D geometry.
- Evaluate whether optimized designs remain effective in 3D.

**Status:** Geometry extrusion utilities exist; a verified 3D flow solver is not yet implemented.

---

## Slide 15 — Verification

### Proposed 3D Verification

- Mesh-independence studies
- Global mass conservation
- Pressure and outlet-flow consistency
- Wall boundary-condition verification
- Comparison with analytical or established benchmark cases

### 2D–3D Comparison

- Outlet flow distribution
- Pressure loss
- Pumping power
- Recirculation and junction flow behavior

**Placeholder:** 3D mesh-convergence results.

---

## Slide 16 — Cont.

### Reconstructed Wall Validation

- Extract solid–fluid interfaces from Brinkman density.
- Compare linear and quadratic interpolation.
- Apply explicit no-slip wall conditions.
- Calculate wall distance for turbulence modeling.
- Compare with a reference explicit-wall simulation.

### Proposed Comparison Metrics

- Outlet flow-rate error
- Pressure-loss error
- Pumping-power error
- Wall-position error
- Mesh sensitivity

**Research question:** Can a density-derived wall reproduce the behavior of an explicitly modeled solid boundary?

---

## Slide 17 — Sensitivities

### 1. Side-Port Flow Sensitivity

\[
G_1=\frac{dQ_{\mathrm{port}}}{d\boldsymbol\gamma}
\]

- Identifies density changes that influence outlet flow distribution.
- Supports flow-routing optimization.

### 2. Pumping-Power Sensitivity

\[
G_2=\frac{dP_{\mathrm{pump}}}{d\boldsymbol\gamma}
\]

- Identifies density changes affecting hydraulic power requirements.
- Supports energy-efficient design.

### 3. Reconstructed Wall-Position Sensitivity

\[
G_3=\frac{\partial y_w}{\partial\boldsymbol\gamma}
\]

- Measures movement of the reconstructed boundary.
- Compares interpolation methods and geometric stability.

### Proposed Methods

- Discrete adjoint differentiation for flow responses
- Direct differentiation for local wall geometry
- Finite-difference and Taylor-test verification

---

## Slide 18 — Optimization

### Primary Design Variables

\[
\boldsymbol\gamma=
[\gamma_1,\gamma_2,\ldots,\gamma_N]^T
\]

- Each density value controls local Brinkman resistance.
- Density changes modify the flow field.
- Adjoint gradients guide optimization updates.

### Proposed Optimization Problem

\[
\max_{\boldsymbol\gamma}\quad Q_{\mathrm{port}}
\]

Subject to:

\[
P_{\mathrm{pump}}\leq P_{\max}
\]

\[
V_{\mathrm{fluid}}\leq V_{\max}
\]

\[
0\leq\gamma_i\leq1
\]

### Optimization Approach

- Gradient-based constrained optimization
- Density filtering and projection
- Iterative flow solution and sensitivity calculation

**Status:** Proposed optimization formulation.

---

## Slide 19 — Cont.

### Additional Potential Controls

**Geometric variables**

- Side-port position
- Side-port width
- Channel dimensions

**Operating parameters**

- Inlet velocity
- Main-outlet pressure
- Side-port pressure

**Numerical/model parameters**

- Maximum Brinkman resistance
- Penalization exponent
- Density filter radius
- Projection sharpness
- Wall-extraction threshold

### Proposed Development Strategy

1. Optimize Brinkman density only.
2. Verify objective and constraint gradients.
3. Investigate additional geometric design variables.
4. Test robustness across operating conditions.
5. Validate candidate designs with explicit-wall models.

---

## Slide 20 — Continuation / Visual

*No established title or text was present in the uploaded slide.*

**Suggested visual content:**

Optimization iteration diagram:

Initial Density → Flow Solver → Objective Evaluation → Adjoint Gradient → Density Update → Convergence Check

**Placeholder:** Optimization history, density evolution, or final geometry once available.

---

## Slide 21 — Summary

### Project Progress

- Established the governing mathematical framework.
- Defined geometry, material, and boundary-condition interfaces.
- Developed numerical benchmark verification tools.
- Verified tangent and adjoint sensitivity calculations on benchmark cases.

### Planned Work

- Complete the fixed-geometry 2D manifold solver.
- Introduce Brinkman density and reconstructed-wall treatment.
- Implement and validate the selected RANS turbulence model.
- Calculate hydraulic and geometric sensitivities.
- Perform topology optimization and 3D validation.

### Main Research Goal

Develop a reliable method for connecting density-based flow optimization with physically meaningful manifold geometry.

---

## Slide 22 — AI Use

### AI-Assisted Development

- Research planning and technical discussion
- Mathematical formulation review
- Programming and debugging assistance
- Documentation and presentation drafting
- Interpretation of numerical verification results

### Verification and Responsibility

- Numerical claims must be supported by reproducible calculations.
- AI-generated code and formulations require independent review.
- Final engineering decisions remain the responsibility of the project author.
- Report AI use in accordance with course requirements.

---

# Presentation Development Notes

- Preserve the existing images on the Background slide.
- Preserve the existing analysis workflow diagram.
- Keep mathematical slides focused on one or two primary equations.
- Replace placeholders with actual simulation figures when available.
- Do not present planned RANS, 3D simulation, or topology optimization as completed.
- Confirm whether the final manifold configuration has two or three outlets; the uploaded presentation specifies three, while a separate proposed side-port study uses two.
- Keep detailed derivations and verification data in supporting documentation rather than overcrowding slides.
