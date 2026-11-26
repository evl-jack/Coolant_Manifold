# Fixed-pipe discrete adjoint verification

The density variables $m$ enter the assembled residual through the Brinkman drag. The mesh, inlet prescription and port definitions are fixed. Eliminate the prescribed velocity degrees of freedom to obtain the free residual $R(z,m)=0$. The reduced state derivative satisfies $R_z z_m=-R_m$. With $R_z^T\lambda_F=\nabla_z F$, the total derivative is

$$
\nabla_m F=\partial_m F-R_m^T\lambda_F.
$$

Both responses depend on the design through the state. Thus $\partial_m F=0$ at fixed $z$. The state Jacobian is the converged Newton tangent restricted to its free rows and columns. Its convection action on $\delta u$ contains $(\delta u\cdot\nabla)u+(u\cdot\nabla)\delta u$. The pressure and continuity couplings follow the solver's paired negative weak-form signs. The adjoint is extended to the prescribed velocity degrees of freedom with zeros.

## Density chain

The editable Helmholtz filter solves $H_{EE}\widetilde\gamma_E=M_E m_E-H_{EF}m_{0,F}$, where $H=M+r_f^2K$. The fixed top neighbors outside the stored grid enter through the eliminated halo across the outlet collar and adjacent solid. The editable forward action is $H_{EE}^{-1}M_E d$. The reverse action is $M_E^T H_{EE}^{-T}s$ for an editable filtered cotangent $s$. The filter keeps the fixed densities binary and can exchange density with the fixed anchors.

For a threshold $t=0.5$ and $\beta>0$, the projected density is

$$
\overline\gamma=\frac{\tanh(\beta t)+\tanh(\beta(\widetilde\gamma-t))}
{\tanh(\beta t)+\tanh(\beta(1-t))},\qquad
\frac{d\overline\gamma}{d\widetilde\gamma}=
\frac{\beta\,\mathrm{sech}^2(\beta(\widetilde\gamma-t))}
{\tanh(\beta t)+\tanh(\beta(1-t))}.
$$

At $\beta=0$ the projection is the identity. The dimensional drag is $\alpha=\alpha_{\max}q(1-\overline\gamma)/(q+\overline\gamma)$ with $d\alpha/d\overline\gamma=-\alpha_{\max}q(q+1)/(q+\overline\gamma)^2$. Every design element uses the constant resistance of its assigned cell at all quadrature points, multiplied by $H_{\mathrm{pipe}}^2/\mu$ with $H_{\mathrm{pipe}}=4\ \mathrm{mm}$. The transpose sums quadrature cotangents back into the corresponding cells with the same scale. The weak drag contributes $\int\alpha^*u\cdot v\,d\Omega$ in the dimensionless equations, so its quadrature cotangent for $F$ is $-w_q u_q\cdot\lambda_{F,q}$ before applying the transfer transpose.

## Responses and finite differences

With outward normals, the reference power convention is $P=-\mu U^2\sum_k\int_{\Gamma_k}(p+Re\,|u|^2/2)(u\cdot n)\,ds$ in nondimensional port coordinates and $U=Re\mu/(\rho H_{\mathrm{pipe}})$. The audit evaluates $f=P/P_{\mathrm{ref}}$ and $\eta=Q_{\mathrm{side}}/Q_{\mathrm{in}}$. The influx $Q_{\mathrm{in}}$ is positive. Its quotient rule is $d\eta=(dQ_{\mathrm{side}}-\eta\,dQ_{\mathrm{in}})/Q_{\mathrm{in}}$. The optimizer uses $g=(\eta_{\mathrm{ref}}-\eta)/0.001$, so $\nabla_m g=-1000\nabla_m\eta$. The same-stage archived immersed original supplies both references.

Each comparison uses a feasible direction with infinity norm one. Interior checks use the centered formula. Near a bound, a direction points into the box and the two forward evaluations at $m+hd$ and $m+2hd$ give a second-order derivative estimate. The code checks the perturbed coordinates directly. Both responses use the same forward states. The independent port evaluator uses five-point edge quadrature. Taylor screening examines $|F(m+hd)-F(m)-h\nabla F\cdot d|$ on consecutive intervals above the recorded screening floor, using the target slope of 1.8. The derivative table retains all step sizes, absolute errors, meaningful relative errors and any failed screens.

## Reproduction

Run from the directory containing `project4_fem`. The first argument must be the directory containing the original and selected `stage_*` snapshots from `Fixed_Pipe_Pilot.zip`. Use one BLAS thread and the same FEM runtime as the pilot.

```bash
python -m project4_fem.project43_topology.verification.fixed_pipe_adjoint_audit PILOT_SNAPSHOTS OUTPUT/adjoint.json
python -m project4_fem.project43_topology.verification.plot_fixed_pipe_adjoint OUTPUT/adjoint.json OUTPUT
```

The JSON records the SHA-256 of each checked snapshot, perturbation directions, the condition number of the editable filter matrix and each Newton tangent comparison. The CSV and the PDF and PNG Taylor plots are generated from that JSON.

