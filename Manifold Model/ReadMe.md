# Fixed-geometry manifold analysis

This implementation lives in `Project 4/Manifold Model`, with a local working
copy retained for development. The existing `Project 4/Solver` and verification
files are separate and untouched. No optimization loop is included.

## Step 1: model and geometry inputs

Implemented:

- `core/geometry/boundary_tags.py`: named inlet, outlet, and wall boundaries.
- `core/geometry/two_d/polygon.py`: explicit counterclockwise polygon outline,
  boundary edge geometry, and outward normals. No dimensions or branch layout
  are supplied by the code. Supports concave outlines with straight boundaries,
  without holes. Ports can span multiple consecutive edges.
- `core/physics/material.py`: supplied constant fluid density and viscosity.
- `core/boundaries/conditions.py`: supplied vector inlet profile, pressure-load
  outlets, and linear-resistance outlets. These describe conditions; finite
  element enforcement is a later step.
- `core/model.py`: connects the inputs for a constant-depth 2D model with one
  inlet and any positive number of outlets. Equal target flow is Qin/n, used
  for reporting and not imposed as an outlet condition.

There are no numerical problem values, example geometries, default fluid
properties, default outlet pressures, or inferred inlet flow rates.
Input validation checks geometry and data consistency, not physical accuracy.

## Established modeling choices

The starting model is steady, incompressible, laminar, constant-property,
depth-averaged flow, with stationary walls and fixed external geometry.
The physical equations and proposed Brinkman interpolation are those in the
existing `Project 4/Outline/Analysis.md`, as revised by the conversation.
No additional body forces, heat transfer, or turbulence assumptions are added.

The fluid's physical density and the Brinkman design field are different
quantities. For the initial hollow interior, the prescribed design field will
be gamma=1, giving zero design drag. Top/bottom plate drag is a separate term
in the established depth-averaged model.

Pressure conditions use a pressure load in the traction condition, not an
additional independently imposed velocity. Linear downstream resistance uses
p_load=p_reference+R*Q, as already specified in the analysis outline. The user
must supply the pressures/resistances. Do not add a zero-mean pressure
constraint automatically when outlet loads already establish the reference.

## Separate dimension-specific folders

- `core/geometry/two_d/` and `three_d/`
- `core/elements/two_d/` and `three_d/`
- `core/rans/two_d/` and `three_d/`

Folders containing only `__init__.py` are reserved locations, not working
solvers. In particular, no RANS or 3D analysis is implemented yet. The 2D
depth-averaged turbulent model needs a stated treatment of the unmeshed depth
before its equations can be implemented. Do not simply attach planar SST
equations to the laminar plate-drag model and claim a verified closure.

## Next coding steps

1. Mesh the explicitly supplied polygon and preserve boundary tags.
2. Implement planar velocity/pressure elements and boundary integration.
3. Add the prescribed Brinkman field and interpolation, with no design updates.
4. Assemble laminar momentum/continuity, inertia, and the established plate drag.
5. Enforce nonzero inlet velocity and outlet traction/resistance conditions.
6. Connect the sparse nonlinear solution and failure reporting.
7. Calculate port flows, pressure-based pumping power, and field exports.
8. Implement fixed-design RANS once its dimensional assumptions are specified.
9. Perform formal benchmark verification and response-gradient checks.

Iterative optimization follows those checks. Later 3D work replaces the
dimension-specific components and removes depth-averaging plate drag while
retaining shared input, solution, and reporting interfaces where applicable.

This staging folder does not yet contain a working flow solver or a runnable
physical case. Numerical values and a particular manifold outline must come
from the problem definition rather than being invented by the implementation.


## Hand references and 3D/sensitivity preparation

- [Symbolic hand reference for all six rungs](docs/verification_hand_reference.md).
- [Researched 3D meshing and gradient plan](docs/three_d_and_gradients.md).
- `core/geometry/three_d/extruded.py` now defines an extrusion of a supplied
  outline, tagged faces, volume, and geometric depth derivatives. The earlier
  description of this geometry package as entirely reserved is superseded.
- `core/sensitivities/local_terms.py` supplies analytic local derivatives of
  the specified Brinkman interpolation and laminar plate drag.
- `core/sensitivities/implicit.py` supplies tangent and adjoint algebra for
  explicitly supplied residual partials; NumPy/SciPy are in requirements.txt.

Only syntax, geometry, and abstract algebra checks have been performed for
these additions. No 3D mesh or 3D flow solve, RANS simulation, complete manifold
response gradient, or iterative optimization is implemented by these additions.
No numerical production inputs have been assigned. The new 3D analysis follows
2D verification and will also need its own verification.
