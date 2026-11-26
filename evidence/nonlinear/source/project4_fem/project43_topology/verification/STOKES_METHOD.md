# Fixed-pipe Stokes state and adjoint

The comparison uses the archived nondimensionalization, viscosity, inlet profile, mesh and traction boundary conditions. At a fixed design, the Taylor–Hood residual is

\[
R(z,m)=A(m)z-b,\qquad A(m)=A_{\mathrm{viscous,pressure}}+A_{\mathrm{Brinkman}}(m).
\]

The original implementation sets the inertia coefficient to zero to omit both convective contributions to the state tangent. The sparse update solves the affine state equation after eliminating prescribed velocity degrees of freedom. The mechanical port objective retains the archived pressure and kinetic-energy flux contributions and the same physical velocity scale. It can be nonlinear in the state even though the governing residual is affine.

For either response \(F\), the prescribed state entries are independent of the raw densities. The existing adjoint solves \(A_{ff}^{T}\lambda=F_{z_f}\) on the free indices and gives

\[
\frac{dF}{dm}=F_m-\lambda^T\frac{\partial A_{\mathrm{Brinkman}}}{\partial m}z.
\]

For fixed ports, the explicit term \(F_m\) is zero. Each quadrature point contributes the negative velocity dot product between the state and the adjoint times its quadrature weight. The existing transpose transfer, rational Brinkman interpolation derivative, projection derivative and anchored Helmholtz transpose map it to the editable raw density variables. The common archived side constraint is \((\eta_{\mathrm{ref,NS}}-\eta)/0.001\).

Run the independent audit from the project root.

```bash
python -B -m project4_fem.project43_topology.verification.fixed_pipe_stokes_audit \
  output/fixed_pipe_pilot_local output/fixed_pipe_stokes_audit/stokes_adjoint_audit.json
```

Optionally include the Stokes endpoint with `--selected path/to/fields.npz`. The audit uses an independently assembled anchored Helmholtz stencil, direct quadrature of each state tangent term, independently integrated port responses, and feasible raw-density perturbations. The original bound design uses second-order inward finite differences. The deterministic gray design uses centered differences. The additional step sizes 0.0005 and 0.0002 probe between the truncation- and roundoff-limited regimes with the original derivative thresholds.
