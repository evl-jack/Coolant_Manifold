# Fixed-pipe Brinkman optimization

Minimize pumping power at the prescribed inlet flow while maintaining the side-outlet flow of the unchanged immersed reference.

| Role | Model |
|---|---|
| Design | 128 dimensionless transition-cell densities bounded by [0, 1] |
| State | P2 velocity and P1 pressure |
| Objective | Mechanical port power normalized by the same-stage reference |
| Constraint | At least the same side-flow fraction as the unchanged immersed manifold |
| Residual | Incompressible nonlinear Navier–Stokes with Brinkman resistance |
| Updates | Discrete adjoint and MMA with β = 0, 2, 8 continuation |

\[
\min_{0\le d\le1}\frac{P(z(d))}{P_0},\qquad
R(z,d)=0,\qquad (\eta_0-\eta)/0.001\le0.
\]

<ol class="loop" aria-label="Optimization loop in execution order">
<li>Cell densities</li>
<li>Anchored Helmholtz filter</li>
<li>Heaviside projection</li>
<li>Rational Brinkman resistance</li>
<li>Nonlinear FEM by Newton</li>
<li>Power and side-flow fraction</li>
<li>Discrete adjoint sensitivities</li>
<li>MMA updates cell densities</li>
</ol>

The updated densities return to the filter. The main pipe and outlet collar retain their original shape throughout the optimization.

\[
\alpha(\gamma)=\alpha_{\max}\frac{q(1-\gamma)}{q+\gamma},
\qquad R_z^T\lambda_J=J_z^T.
\]

[Official Project 4 assignment](https://designinformaticslab.github.io/DesignOptimization2025/project4.html)
