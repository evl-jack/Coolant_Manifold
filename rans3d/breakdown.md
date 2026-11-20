# In-house 3D coolant manifold solver

## 1. Scope and lineage

This solver advances the three mean velocity components, pressure, turbulent kinetic energy and specific dissipation rate. It solves the steady incompressible equations on a cell-centered finite-volume mesh. The laminar mode removes the turbulence transport equations and sets eddy viscosity to zero. The 3D code runs with Python, NumPy and SciPy. It can read the repaired Sol mesh in OpenFOAM’s binary format. File-format compatibility lets us use the same geometry for comparison while running our own matrix assembly, pressure correction and turbulence transport.

The 2D code introduced first-order upwind transport, least-squares gradients, Rhie–Chow pressure interpolation and SIMPLE-style pressure correction on a Cartesian grid. It also uses Brinkman penalization to explore changes to the passage shape. The 3D solver uses the selected shape’s body-fitted wall patches. It solves directly for $u_x,u_y,u_z$ and has one volume per cell and one area vector per face. The repaired 3D mesh comes from the shape inferred in the 2D design search. [2D solver](../rans_solver/solve.py) · [2D Brinkman design notes](../rans_solver/BRINKMAN.md) · [Gradient and mesh repair notes](../GRADIENTS_2D_AND_3D.md)

| File | Owns |
|---|---|
| [mesh.py](mesh.py) | Sol binary mesh import, polygon and cell geometry, connectivity, Cartesian test duct |
| [solver.py](solver.py) | Face interpolation, gradients, matrix assembly, pressure correction, SST-style equations, convergence tests |
| [__main__.py](__main__.py) | Command line, Poiseuille reference solution, restart input, histories and field output |

## 2. Equations and units

For constant-density water, let $\boldsymbol u$ be mean velocity in m/s, $\nu=\mu/\rho$ molecular kinematic viscosity in m²/s, and $\nu_t$ turbulent kinematic viscosity in m²/s. The solver uses modified kinematic pressure $p^*=p_{\mathrm{static}}/\rho+2k/3$. This absorbs the isotropic Reynolds stress under the Boussinesq eddy-viscosity assumption.

$$
\nabla\cdot\boldsymbol u=0,\qquad
\nabla\cdot(\boldsymbol u\otimes\boldsymbol u)
=-\nabla p^*+\nabla\cdot[(\nu+\nu_t)(\nabla\boldsymbol u+\nabla\boldsymbol u^\mathsf{T})]
$$

The implementation handles the scalar Laplacian of each velocity component implicitly and the transposed-gradient contribution explicitly on internal faces. Its explicit transposed-stress flux is zero on boundary faces, which is a boundary approximation. Pressure saved to disk is static gauge pressure in Pa, obtained through $p_{\mathrm{static}}=\rho(p^*-2k/3)$.

With $S=\sqrt{2S_{ij}S_{ij}}$, where $S_{ij}$ is the symmetric mean strain tensor, the SST-style closure limits eddy viscosity and production through

$$
\nu_t=\frac{a_1 k}{\max(a_1\omega,SF_2)},\qquad
P_k=\min(\nu_tS^2,10\beta^* k\omega),\qquad
a_1=0.31,\ \beta^*=0.09.
$$

The transported variables are $k$ in m²/s² and $\omega$ in s⁻¹. In kinematic form, their steady equations are

$$
\begin{aligned}
\nabla\cdot(\boldsymbol u k)
 &=\nabla\cdot[(\nu+\sigma_k\nu_t)\nabla k]+P_k-\beta^* k\omega,\\
\nabla\cdot(\boldsymbol u\omega)
 &=\nabla\cdot[(\nu+\sigma_\omega\nu_t)\nabla\omega]
   +\gamma P_k/\nu_t-\beta\omega^2
   +2(1-F_1)\frac{\sigma_{\omega2}}{\omega}\nabla k\cdot\nabla\omega.
\end{aligned}
$$

The code computes $F_1$ and $F_2$ from $k,\omega,\nu$ and the estimated distance to the wall. $F_1$ blends near-wall and outer coefficients. $F_2$ enters the viscosity limiter. Destruction terms are implicit with lagged coefficients. Positive cross-diffusion is a source and negative cross-diffusion is treated as a sink. Lower bounds of $10^{-12}$ for $k$ and $10^{-8}$ for $\omega$ guard positive denominators. The turbulence model is an initial implementation based on SST-2003. The wall treatment and distance approximation require validation for the manifold. [Turbulence Modeling Resource SST equations](https://tmbwg.github.io/turbmodels/sst.html) · [SST implementation](solver.py)

## 3. Mesh geometry and boundary conditions

The imported mesh is an undecomposed, single-region binary `polyMesh` with compact face storage, little-endian 32-bit labels and 64-bit coordinates. Internal faces come first. `owner` and `neighbour` identify their two cells. Boundary faces have only an owner and a patch name. Coordinates, face area vectors and cell volumes are expressed in SI units.

1. Split each polygon into triangles about its mean vertex position. Sum triangle area vectors to obtain the oriented face area vector. Use projected triangle weights for the face centroid.
2. Use signed tetrahedra from each face to a cell reference point to accumulate cell volume and volume centroid. Each internal face enters both neighboring cells with opposite signs.
3. Check positive volumes, boundary addressing, one connected cell graph and positive projected owner-to-neighbor or owner-to-boundary distances. The solver expects the mesh to have passed the full mesh quality check separately.
4. Precompute face-normal distances, interpolation weights and least-squares stencils. Near-wall SST distance is approximated with the nearest wall face center. Wall-adjacent cells also use their projected distance to their wall faces. [Mesh reader and duct generator](mesh.py) · [Mesh quality workflow](../sol_rans_3d/README.md)

| Patch | Mean velocity | Pressure | Turbulence |
|---|---|---|---|
| `inlet` | Fixed vector. The manifold default is 2 m/s along x. | Inlet flux is prescribed by velocity. Pressure gradient reconstruction uses interior cells. | Fixed $k$ and $\omega$ from inlet speed, intensity and length scale. |
| `main_outlet`, `side_outlet` | Owner-cell velocity for predicting normal flux. The pressure correction updates the flux. | Zero static gauge pressure. The known modified pressure uses $2k/3$ from the owner. | Outflow uses owner values. Reverse convection uses the inlet turbulence values. |
| `walls` | No slip with exactly zero normal volumetric flux. | No correction to wall flux. | $k=0$, with $\omega=6\nu/(\beta_1 d_n^2)$ as a face-distance viscous estimate. |

The default water values are $\rho=998.2$ kg/m³ and $\mu=1.002\times10^{-3}$ Pa·s. With the default manifold inlet, $I=0.05$ and $\ell=0.56$ mm give $k_\mathrm{in}=1.5(U I)^2=0.015$ m²/s² and $\omega_\mathrm{in}=\sqrt{k_\mathrm{in}}/(0.09^{1/4}\ell)\approx399.3$ s⁻¹. Intensity and length scale can be varied independently from the command line. They are modeling assumptions inherited from the 2D baseline. The transport boundary takes zero velocity for the incoming convective term on an outlet face and reports the amount of reverse flow for diagnosis. Strong outlet backflow warrants a check of the modeled outlet location. [Settings and boundary implementation](solver.py)

## 4. What happens in a solver iteration

1. **Reconstruct gradients.** For each cell, fit the neighboring value differences to their position differences using inverse squared distance as the weight. Include known boundary values where applicable. The resulting 3 × 3 geometric inverse is cached for each stencil.
2. **Update lagged transport coefficients.** In SST mode, evaluate strain, blending, eddy viscosity, limited production and cross-diffusion. Add $\nu_t$ to the momentum diffusion coefficient and use the blended diffusivities in the turbulence equations.
3. **Build momentum systems.** A face flux uses the upstream cell’s value in the implicit convection matrix. Conductance is $\Gamma_f A_f/d_n$, where $d_n$ is the center-to-center or center-to-boundary distance projected on the face normal. Internal diffusion adds an explicit nonorthogonal correction. The three momentum components share one matrix. A local pseudo-time term $V(q^{n+1}-q^n)/\Delta t$ stabilizes the iteration, with $\Delta t=\mathrm{CFL}\,V^{1/3}/U_\mathrm{ref}$. At a fixed point this term vanishes.
4. **Predict and correct flux.** The pressure interpolation subtracts the part of the pressure jump already represented by the interpolated gradient. This Rhie–Chow-style correction couples face velocity to neighboring pressures on the collocated mesh. Let $F_f^*$ be the predicted flux and $D_f=(V/a_P)_f A_f/d_n$. The pressure-correction equation adjusts each flux by $D_f(p'_P-p'_N)$, with $p'=0$ on the pressure boundaries. The same corrected face flux is added to the owner and subtracted from the neighbor. Cell velocity receives the pressure-gradient correction, and cell pressure uses relaxation factor 0.3.
5. **Update $k$ and $\omega$.** Assemble their advection-diffusion systems using the corrected flux and the lagged SST coefficients. The new values feed the next momentum iteration.
6. **Record diagnostics.** Save the largest normalized change in the solved fields, local and global continuity, inlet pressure, outlet flows, reverse flows and the running variation of the integral quantities. [Finite-volume implementation](solver.py)

Convection order 1 uses implicit upwind throughout. Order 2 keeps that matrix and adds a deferred linear-upwind correction. The extrapolated internal face value is clipped to the range of the adjacent cells. This gives the higher-order reconstruction where the limiter allows it. The limiter can locally reduce the effective order. Duct runs use sparse LU. Meshes above 8,000 cells use diagonally preconditioned BiCGSTAB for transport and conjugate gradient for the symmetric pressure correction. The Krylov tolerance is $10^{-7}$. The pseudo-time term advances an iterative steady calculation. Time-accurate URANS would require physical time integration and separate verification.

## 5. Convergence and outputs

The default steady solve requires all of these conditions:

- The maximum normalized field change is below $10^{-5}$. Velocity and modified pressure changes are scaled with $U_\mathrm{ref}$ and $U_\mathrm{ref}^2$, and turbulence changes use the inlet or current field scales.
- The sum of absolute cell flux defects and the absolute net boundary flux are each below $10^{-5}$ of inlet flow. Each is normalized separately.
- After 20 samples, each outlet flow’s range over the window is at most 0.1% of inlet flow. The inlet pressure drop range is at most 0.1% of its absolute mean.

Iteration count limits total work. `summary.json` reports whether these conditions were met, the final metrics and numerical settings. `history.csv` contains per-iteration diagnostics, with positive inlet magnitude and signed outward outlet flow. `fields.npz` contains centers in m, velocity in m/s, static gauge pressure in Pa, $k,\omega,\nu_t$, and oriented face fluxes in m³/s. Summary patch flows are signed outward, so inlet flow is negative. Reverse flow is the positive magnitude of inward flux on individual outlet faces. The reported pressure drop uses face-area-weighted inlet-owner pressure relative to the zero-pressure outlet boundary. The duct verification also fits the axial pressure gradient in an interior section to avoid depending on inlet cell sampling.

The residual used by a pressure linear solve measures the discrete correction equation. Small continuity error can therefore occur early in an unconverged calculation. The field-change and integral monitors add evidence of a steady fixed point. Grid independence, turbulence model validation and the manifold design target are separate checks. The design target is $Q_\mathrm{side}/Q_\mathrm{main}=0.75$. [Monitoring and output implementation](__main__.py)

## 6. Verification performed

The analytical test uses a 20 mm long square duct of side $H=3.5$ mm. For a fully developed laminar duct, the velocity satisfies $-\nu(\partial_{yy}+\partial_{zz})u_x=G$ with zero velocity on all four walls. A double sine series gives the profile and a reference pressure gradient of **116.372 Pa/m** at a continuous cross-sectional mean speed of 0.05 m/s. The exact inlet profile is sampled at the inlet faces. Discrete midpoint quadrature makes the prescribed net flux depend slightly on grid resolution. We compare velocity and fitted pressure gradient in the middle half of the duct.

| Case | Cells | Iterations | Relative velocity L2 error | Pressure gradient error |
|---|---:|---:|---:|---:|
| Laminar, 12 × 12 × 12 | 1,728 | 55 | 0.819% | 3.80% |
| Laminar, 12 × 24 × 24 | 6,912 | 70 | 0.206% | 1.50% |
| Laminar, 12 × 12 × 12, order-2 restart | 1,728 | 21 after restart | 0.798% | 3.91% |

The transverse refinement lowered both measured errors. These two resolutions check the trend. A formal order-of-accuracy study would use more grids and consistent refinement in all directions. The 1,728-cell SST-style duct run also converged with the limited linear-upwind option in **108 iterations**, with approximately **2.63 kPa** inlet-to-outlet pressure drop. This verified that both turbulence transport equations and the higher-order path ran to a steady numerical state on that mesh. The SST duct pressure drop has no independent benchmark here.

The solver also imported the **306,218-cell repaired manifold** and completed five pressure-corrected SST iterations with normalized net mass imbalance below $2\times10^{-8}$ throughout. This exercised the polyhedral import, addressing, patch fluxes and large-mesh sparse solver path. The initial pressure estimate changed from +358 kPa to −201 kPa between the first two iterations and still fluctuated at iteration five. Inlet pressure and outlet split are startup diagnostics until the manifold converges. The solver currently uses one Python process. Distributed MPI and pressure multigrid are future performance work. [Reference profile and reporting](__main__.py) · [Validated source mesh](../MANIFOLD_HANDOFF.md#6-hex-mesh-repair-and-sol-restart)

## 7. Reproduce and continue a run

From the repository root in the existing Windows virtual environment:

```powershell
.\.venv\Scripts\python.exe -m rans3d --laminar --output runs/duct12
.\.venv\Scripts\python.exe -m rans3d --laminar --ny 24 --output runs/duct24
.\.venv\Scripts\python.exe -m rans3d --order 2 --output runs/duct-sst
.\.venv\Scripts\python.exe -m rans3d --laminar --order 2 --restart runs/duct12/fields.npz --output runs/duct12-order2
```

For a new manifold input, give the path to the undecomposed repaired mesh. Start with a short run to check mesh connectivity and finite fields. Retain the first run’s `fields.npz` if the higher-order run will use it as a checkpoint.

```powershell
.\.venv\Scripts\python.exe -m rans3d --mesh path/to/polyMesh --steps 5 --output runs/manifold-smoke
```

The checkpoint needs the same centers, cell ordering, face ordering and turbulence mode. It stores cell fields and face fluxes, so first-order fields can initialize a separate order-2 run. The current in-house CLI leaves the switch decision to the caller. Compare changes and monitor histories before switching. The separate OpenFOAM case retains its automatic staged test.

On Sol, use a Python environment with NumPy and a SciPy version that supports `rtol` in its iterative solvers. Run `python -m rans3d` from a checkout containing both `rans3d/` and the shared binary reader in `sol_rans_3d/repair_hex_mesh.py`. Run it on an allocated compute node and pass the path to `repaired-coarse/constant/polyMesh`. The current implementation uses a single process. The fine mesh needs a separate memory and runtime check before a full solve.

