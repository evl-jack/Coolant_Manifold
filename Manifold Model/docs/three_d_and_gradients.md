# 3D implementation and sensitivity plan

Status: researched design and initial geometry/derivative building blocks.
No 3D flow solution, mesh-quality result, RANS result, or complete PDE gradient
is claimed. No optimization loop is introduced. The first 3D analysis is
scheduled after the 2D model and its derivatives pass their verification.

## Recommended path

My engineering recommendation is to first extrude the supplied 2D outline
through its actual constant depth. It retains the same port names and supports
any number of outlets without inventing a different manifold. It resolves
top/bottom walls and allows nonzero out-of-plane velocity. This is an initial
geometry family, not a restriction on future imported CAD, round pipes,
fillets, or manifolds with ports on different planes.

Gmsh supports extrusion, CAD operations, mesh-size fields, and physical groups.
Use named surface groups for every port and wall, plus a named fluid volume;
rebuild semantic groups after CAD Boolean operations rather than assuming
numeric entity IDs stay unchanged. My initial mesh recommendation is a
conforming tetrahedral mesh with local refinement. A later swept prism/hex
mesh can reduce cost in suitable geometries, but requires its own compatible
element implementation. The 2D BoundaryLayer size field is not a general
3D layer generator; 3D layers require an appropriate extrusion/block strategy.
[Gmsh reference manual](https://gmsh.info/doc/texinfo/gmsh.html)

The 3D laminar model uses physical volume and surface integrals:

$$
\nabla\cdot\mathbf u=0,\qquad
\rho(\mathbf u\cdot\nabla)\mathbf u
-\nabla\cdot[2\mu\varepsilon(\mathbf u)]
+\nabla p+\alpha(\gamma)\mathbf u=0.
$$

There is no 12*mu/b^2 drag in this wall-resolved 3D equation: keeping it
would count the top/bottom wall resistance twice. Initially the hollow
interior has gamma=1, hence zero design drag. Keep support for prescribed
Brinkman fields as a separate model term, with no design updates.

## Discretization and numerical efficiency

Keep the initial finite element family consistent with the existing work:
P2 vector velocity and P1 pressure. A tetrahedron has ten local scalar P2
velocity nodes and four P1 pressure nodes: 34 local flow unknowns before
sharing nodes and applying conditions. Check stability on the selected mesh
family; element names alone do not establish accuracy or robustness.

Use curved geometry representations if later circular/filleted CAD is to
be resolved accurately; merely increasing solution order does not correct a
coarse faceted wall. Mesh generation and element import must agree on local
node ordering and face orientation.

Precompute basis data at quadrature points. Cache element mapping data for
the fixed mesh. Reuse sparse structure during nonlinear iterations. Expose
residual evaluation separately from Jacobian assembly so line searches do not
build an unused Jacobian. Never convert the full sparse matrix to dense just
to calculate a diagnostic norm.

For large problems I recommend a replaceable PETSc solver backend: GMRES for
the nonsymmetric linearized system and a velocity/pressure field split with
a problem-appropriate Schur approximation. This is a proposed starting point,
not a claim of mesh-independent convergence. Tolerance, Reynolds number,
Brinkman contrast, and preconditioner quality still matter.
[PETSc field-split documentation](https://petsc.org/release/manualpages/PC/PCFIELDSPLIT/),
[GMRES documentation](https://petsc.org/release/manualpages/KSP/KSPGMRES/)

The official DOLFINx Stokes example demonstrates Taylor-Hood spaces and block
solvers. UFL also supports symbolic differentiation of variational forms.
These are useful references or a future backend, not evidence that legacy
dolfin-adjoint automatically supports this custom solver or a particular
DOLFINx release. Check compatibility before adopting a stack, and do not
silently replace required student-written finite element code.
[Taylor-Hood demo](https://docs.fenicsproject.org/dolfinx/main/python/demos/demo_stokes.html),
[UFL form language](https://docs.fenicsproject.org/ufl/main/manual/form_language.html)

## Mesh plan

1. Accept actual dimensions and labeled geometry; no production mesh can be
   selected without them. Create the 3D fluid volume from the planar outline
   initially, or later from explicitly supplied CAD.
2. Preserve named inlet/outlet/wall groups into the mesh import. Require every
   exterior facet to have an intended boundary role and every element to have
   a supported cell type. Confirm ports have nonzero area and fluid is connected.
3. Refine branch entrances, junctions, narrow passages, and physical walls;
   avoid spending the same resolution everywhere. Feature sizes and requested
   accuracy, not an arbitrary element count, set the mesh controls.
4. Measure positive element Jacobians, shape quality, local size transitions,
   port areas, and normal orientation. Reject invalid cells. For curved
   high-order cells, check Jacobians inside cells, not only at vertices.
5. Check the number and distribution of elements through depth. For a thin
   manifold, isotropic tetrahedra may be costly. A swept anisotropic mesh is
   a later candidate if its quality, conditioning, and elements are supported.
6. Generate a mesh family and assess convergence of each outlet flow and
   pressure-based power. Good-looking cells alone do not establish solution
   accuracy. Record memory and solve times with errors.

These are my proposed quality requirements, not completed mesh tests or
universal numerical thresholds. Mesh-size choices remain required inputs.

## RANS scope

For 3D fixed-design assessment use a precisely named SST variant, initially
the planned SST-2003. Its eddy viscosity and blending depend on the state and
wall distance, and the limited production enters both transport equations.
Inlet k/omega and the wall treatment must be specified. Implement the selected
variant faithfully; choosing SST-2003m or smoothing a limiter is a documented
change, not an interchangeable implementation detail.
[Turbulence Modeling Resource: SST](https://tmbwg.github.io/turbmodels/sst.html)

For a wall-resolved strategy, set the first off-wall spacing from
y_plus=rho*u_tau*y/mu and a stated target; revise it from the computed shear.
There is no defensible universal dimensional first-layer height without
flow and geometry inputs. Layer count, growth, and thickness must resolve
the flow and fit the local geometry. Check sensitivity of the port results.

The fixed hollow geometry has explicit walls. A later Brinkman field does
not automatically provide the wall location or turbulence boundary conditions
of an immersed solid. Do not invent SST wall distance or damping rules for
porous solids. Resolve a sharp extracted geometry for RANS assessment, or
adopt and verify a documented turbulence/porosity treatment separately.

If a physically unsteady junction flow does not admit a useful steady
solution, assess a transient model rather than interpreting failed steady
iterations as reliable physics. That extension and its time-averaged
sensitivities are outside the current stationary implementation.

## What 'a gradient for every parameter' means

For each response J, compute dJ/dm_j for every declared independent continuous
input m_j, including fixed inputs if sensitivity is requested. This is not
the spatial gradient grad(u), nor is it one derivative independent of response.
The full output is a response-by-parameter sensitivity map. Outlet count,
mesh connectivity, model dimension, and turbulence-model selection are
discrete choices and have no ordinary infinitesimal gradient. Compare those
as separate models instead. Numerical tolerances are convergence controls,
not physical design variables.

Derived quantities follow the chain rule. For example, if viscosity is
derived from temperature, differentiate through that property relation if
it is part of the actual model. Do not make both quantities independent
controls without stating the intended experiment. No temperature dependence
is currently supplied and none is invented here.

| Parameter family | Derivative pathway required |
|---|---|
| Fluid density | Inertia; also turbulence and any explicitly density-dependent inlet/measurement definitions |
| Dynamic viscosity | Viscous stress; 2D plate drag; SST molecular diffusion and wall treatment |
| Inlet-profile coefficients | Nonzero Dirichlet lifting, Qin, and any dependent turbulence inlet values |
| Outlet pressure loads/reference pressure | Traction residual and explicit pressure-based response terms |
| Each downstream resistance | Nonlocal outlet load p_ref+R_i*Q_i and its state Jacobian |
| Manifold depth | In 2D: drag, integration measures, inlet/outlet flow and loads; in 3D: geometry and mesh motion |
| Header/branch dimensions, port positions, fillet radii | CAD parameterization, mesh coordinates, mapped derivatives, integration measures, normals, and moving boundary data |
| Prescribed gamma coefficients | Brinkman drag and the field representation; initial field can still be fixed |
| Alpha_max and q_alpha | Drag interpolation; include dependence of alpha_max on other inputs if defined that way |
| Filter/projection controls, if later introduced | Filter solve and projection chain; absent controls must not be reported as implemented |
| Inlet k and omega, if RANS is active | Turbulence boundary conditions and all coupled mean-flow/turbulence equations |
| Turbulence coefficients, if deliberately varied | Full selected closure, including blends and active limiter branches; normally hold published calibration fixed |
| Reference pressure used only for reporting | Explicit derivative in the chosen response; distinguish from a physical outlet load |

For pressure-based power, state the pressure reference and system boundary:

$$
P_{hyd}=Q_{in}(p_{in,fw}-p_{ref}).
$$

This follows the existing project definition and includes downstream pressure
loss when p_ref is the downstream reference. For manifold-only pressure work,
use inlet pressure work minus the sum of outlet pressure works instead.
Electrical pumping power additionally requires a specified pump efficiency.
Do not treat these quantities as identical when outlet pressures differ.

For each outlet define smooth signed performance deviations

$$
g_i=\frac{nQ_i}{Q_{in}}-1,\qquad
\frac{dg_i}{dm}=\frac{n}{Q_{in}}\frac{dQ_i}{dm}
-\frac{nQ_i}{Q_{in}^2}\frac{dQ_{in}}{dm}.
$$

Retain Q_in derivatives when inlet velocity or geometry changes. Report
both signed deviations rather than differentiating a nonsmooth maximum of
absolute deviations. Future upper/lower constraints can use +/-g_i-tolerance.

## Residual-based derivatives

Let z collect all free flow degrees of freedom, and k/omega if active.
After boundary enforcement, R(z,m)=0. For geometry controls, write mesh
coordinates X=X(m) explicitly, so

$$
R_m=\left.\frac{\partial R}{\partial m}\right|_{z,X}
+R_X X_m.
$$

The following forms use the total explicit parameter partial, including
geometry effects and boundary lifting, while holding the free state fixed:

$$
A=R_z,\qquad Az_{m_j}=-R_{m_j},
\qquad A^T\lambda=J_z^T,
\qquad dJ/dm_j=J_{m_j}-\lambda^TR_{m_j}.
$$

Tangent solves suit a small number of parameters and many responses. Adjoint
solves suit many parameters and relatively few responses. Evaluate field
derivatives locally or as vector-Jacobian products rather than materializing
a dense global state-by-field Jacobian.

For nonzero prescribed inlet velocity, include derivatives of the eliminated
values and lifting. For resistance outlets include the derivative of the
integrated Q_i in the state Jacobian, not only an explicit derivative of R_i.
For shape changes include both the surface measure and normal variation in
Q_i. Define whether the inlet profile stays attached to reference port
coordinates or fixed spatial coordinates when a port moves.

Use fixed-connectivity mesh deformation for shape derivative checks. A fresh
remesh changes the discrete problem and can contaminate finite differences.
Do not differentiate mesh generation decisions as if they were smooth CAD
parameters. Remesh when needed for quality, then recompute and recheck the
derivative on the new mesh.

For RANS, a full derivative includes k, omega, mu_t, blending, wall distance
when geometry changes, and the selected wall condition. Frozen-turbulence
derivatives are an approximation and must be labeled as such. min/max
limiters and nearest-wall switches are piecewise differentiable. Check smooth
branches; at a switch use directional/generalized reasoning or a documented
smoothing, not a claim of a unique classical gradient everywhere.

## Analytic local derivatives already available for implementation

For alpha=alpha_max*(1-gamma)/(1+q*gamma),

$$
\alpha_\gamma=-\frac{\alpha_{max}(1+q)}{(1+q\gamma)^2},\quad
\alpha_{\alpha_{max}}=\frac{1-\gamma}{1+q\gamma},\quad
\alpha_q=-\frac{\alpha_{max}\gamma(1-\gamma)}{(1+q\gamma)^2}.
$$

For c_plate=12*mu/b^2,

$$
(c_{plate})_\mu=12/b^2,\qquad (c_{plate})_b=-24\mu/b^3.
$$

For V=A*b, V_b=A at fixed planar shape. If the extrusion bottom stays fixed,
top vertex z coordinates have derivative one with respect to b; bottom
vertices have derivative zero. These local expressions do not by themselves
constitute total flow or pumping-power gradients.

## Deliverables and order

1. Complete fixed-geometry 2D forward analysis.
2. Derive and execute the six-rung verification and applicable 2D gradients.
3. Build the extruded CAD/mesh bridge with persistent semantic boundary names.
4. Implement tetrahedral integration, 3D boundary conditions and laminar flow.
5. Verify the new 3D mathematical components and responses independently.
   A passed 2D model does not automatically verify the new 3D implementation.
6. Implement full 3D parameter/shape derivatives and recheck Taylor behavior.
7. Add the selected RANS closure and coupled gradients with its own benchmark
   and near-wall sensitivity assessment.

Keep a registry recording each parameter's units, independence/dependencies,
derivative implementation, affected residual terms, responses, and check
status. Missing derivatives should fail explicitly, never silently return
zero. Optimization remains disabled until the necessary checks pass.
