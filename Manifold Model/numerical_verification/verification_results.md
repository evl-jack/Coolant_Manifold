# Manifold model verification results

This report records the independent analytical references for the six verification rungs. The decimal values are evaluations of closed-form expressions. The solver measurement and pass/fail fields are pending.

## Shared model and inputs

The flow is steady, incompressible, and depth-averaged. Use $\rho=1000\ \mathrm{kg/m^3}$, $\mu=0.001\ \mathrm{Pa\,s}$, and $b=0.005\ \mathrm m$. Let

$$
K_p=\frac{12\mu}{b^2},\qquad
\alpha=\alpha_{\max}\frac{1-\gamma}{1+q_\alpha\gamma},\qquad K=K_p+\alpha.
$$

Here $\alpha_{\max}=480\ \mathrm{Pa\,s/m^2}$ and $q_\alpha=2$. Plate drag is active unless a case explicitly sets total drag to zero. Square benchmarks use uniform $\gamma=1$, so $\alpha=0$ and $K=480\ \mathrm{Pa\,s/m^2}$.

## 1. Exact strip

The strip has $L=0.10\ \mathrm m$, $H=0.020\ \mathrm m$, no slip at $y=0,H$, periodic ends, zero-mean pressure, and $f=(G,0)$ with $G=1.0\ \mathrm{N/m^3}$. The fully developed solution has $U=(u(y),0)$ and $p=0$. The equivalent driving pressure drop is $GL=0.10\ \mathrm{Pa}$.

Set $s=H/2$, $\lambda=\sqrt{K/\mu}$, and $z=\lambda s$. The equation $\mu u''-Ku=-G$ and the wall values give

$$
u(y)=\frac{G}{K}\left[1-\frac{\cosh(\lambda(y-s))}{\cosh z}\right],\quad
Q=\frac{bGH}{K}\left(1-\frac{\tanh z}{z}\right),\quad
\tau_w=\frac{G}{\lambda}\tanh z,\quad \mathcal P=GLQ.
$$

The signed streamwise traction exerted on the fluid is $-\tau_w$ at each wall.

| Analytical quantity | Hollow, $\gamma=1$ | Active Brinkman drag, $\gamma=0.5$ |
|---|---:|---:|
| $K_p$ (Pa s/m²) | 480 | 480 |
| $\alpha$ (Pa s/m²) | 0 | 120 |
| $K$ (Pa s/m²) | 480 | 600 |
| $z$ | 6.928203 | 7.745967 |
| Centerline velocity (m/s) | $2.079251\times10^{-3}$ | $1.665225\times10^{-3}$ |
| Volumetric flow (m³/s) | $1.782631\times10^{-7}$ | $1.451501\times10^{-7}$ |
| Side-wall shear magnitude (Pa) | $1.443373\times10^{-3}$ | $1.290994\times10^{-3}$ |
| Driving power (W) | $1.782631\times10^{-8}$ | $1.451501\times10^{-8}$ |

With both drag terms disabled in a separate benchmark, $u_0=Gy(H-y)/(2\mu)$. The centerline speed is $0.0500\ \mathrm{m/s}$, $Q_0=3.333333\times10^{-6}\ \mathrm{m^3/s}$, $\tau_{w,0}=0.0100\ \mathrm{Pa}$, and $\mathcal P_0=3.333333\times10^{-7}\ \mathrm W$.

## 2. Stokes mesh convergence

Use $[0,S]^2$, $S=0.10\ \mathrm m$, $X=x/S$, $Y=y/S$, $V=0.001\ \mathrm{m/s}$, and $p_*=0.001\ \mathrm{Pa}$. Use zero-mean pressure and prescribe $(u,v)=(0,0)$ on each side. Define

$$
F(t)=t^2(1-t)^2,\quad F'=2t-6t^2+4t^3,\quad
F''=2-12t+12t^2,\quad F'''=-12+24t.
$$

$$
\begin{aligned}
A &= (F(X)F'(Y),-F'(X)F(Y)),\\
B &= (F''(X)F'(Y)+F(X)F'''(Y),-F'''(X)F(Y)-F'(X)F''(Y)),\\
U_2 &= VA,\qquad p_2=p_*(X-\tfrac12)(Y-\tfrac12),\\
f_2 &= 0.01(Y-\tfrac12,X-\tfrac12)-0.0001 B+0.48 A .
\end{aligned}
$$

The divergence terms cancel and $F=F'=0$ at the endpoints. The force is in N/m³ and equals $\nabla p_2-\mu\Delta U_2+KU_2$. Disable convection for this rung.

Use P2 velocity and P1 pressure on triangular meshes with 8, 16, and 32 subdivisions per side if supported. Record actual maximum element diameters and counts. Under the stated smooth Stokes setup with stable meshes and adequate integration, the asymptotic orders for velocity $L^2$, velocity gradient $L^2$, and gauge-aligned pressure $L^2$ are 3, 2, and 2, respectively. Calculate each observed order as $\log(E_c/E_f)/\log(h_c/h_f)$.

## 3. Global balance

The periodic face flows are $-Q$ and $+Q$. Across the strip, integrating the momentum equation gives $GbHL=2bL\tau_w+K_pLQ+\alpha LQ$. Multiplying that equation by $u$ and integrating by parts gives $\mathcal P=bL(\mu I_1+K_p I_0+\alpha I_0)$, where $I_0=\int_0^H u^2\,dy$ and $I_1=\int_0^H(u')^2\,dy$. For direct calculation,

$$
I_0=(G/K)^2[H-3\tanh z/\lambda+s\,\operatorname{sech}^2 z],\qquad
I_1=(G/K)^2\lambda^2[\tanh z/\lambda-s\,\operatorname{sech}^2 z].
$$

| Analytical balance term | Hollow | Active Brinkman drag |
|---|---:|---:|
| Driving force (N) | $1.000000\times10^{-5}$ | $1.000000\times10^{-5}$ |
| Side-wall resistance, both walls (N) | $1.443373\times10^{-6}$ | $1.290994\times10^{-6}$ |
| Plate resistance (N) | $8.556627\times10^{-6}$ | $6.967205\times10^{-6}$ |
| Brinkman resistance (N) | 0 | $1.741801\times10^{-6}$ |
| In-plane viscous dissipation (W) | $1.503473\times10^{-9}$ | $1.075822\times10^{-9}$ |
| Plate dissipation (W) | $1.632283\times10^{-8}$ | $1.075135\times10^{-8}$ |
| Brinkman dissipation (W) | 0 | $2.687838\times10^{-9}$ |

For an open manifold, $\sum_j Q_j=Q_{\rm in}$ with positive inward inlet flow and signed outward outlet flows. An open power check also includes inlet reaction traction work, outlet traction work, and net outward kinetic-energy flux. For full uniform outlet traction $t=-p_j n$, outlet work is $-p_jQ_j$. Numerical open-domain power targets require the outlet traction convention and any resistance law.

## 4. Rotation

Rotate the Rung 2 mesh and all problem data $90^\circ$ counterclockwise about the origin using $R(a,b)=(-b,a)$. The rotated domain is $[-S,0]\times[0,S]$. At $(x',y')$, use

$$
U_4=R U_2(y',-x'),\qquad p_4=p_2(y',-x'),\qquad f_4=R f_2(y',-x'),\qquad n'=Rn.
$$

The prescribed velocity remains zero on all four mapped sides. The pressure reference stays zero mean, and the linear Stokes operator and drag coefficients carry over. Compare $R^\mathsf T U_{4,h}(Rx)$ and $p_{4,h}(Rx)$ to their Rung 2 counterparts. Corresponding scalar error norms agree under the rotation.

## 5. Nonlinear benchmark

Use the same square, scales, material coefficients, and pressure reference. Retain convection. The exact fields and the body force (in N/m³) are

$$
\begin{aligned}
U_5&=V(X^2,-2XY), &p_5&=p_*(X-Y),\\
f_5&=(0.02X^3+0.0098+0.48X^2,\ 0.02X^2Y-0.01-0.96XY).
\end{aligned}
$$

Prescribe the exact trace on all four sides. In order $x=0$, $x=S$, $y=0$, $y=S$, the values are $(0,0)$, $V(1,-2Y)$, $V(X^2,0)$, and $V(X^2,-2X)$.

The convective acceleration is $(V^2/S)(2X^3,2X^2Y)$. The velocity and pressure are representable in P2/P1 on affine triangles with exact trace data and sufficient quadrature. The directional Newton derivative is $\rho[(W\cdot\nabla)U+(U\cdot\nabla)W]$. Use admissible homogeneous velocity increments for a Jacobian-action check. Near the root, a nonsingular Jacobian and accurately solved full Newton steps give local quadratic convergence.

## 6. Strip response sensitivities

Use the $\gamma=0.5$ strip. Define $C=(bH/K)(1-\tanh z/z)$ and $T=\tanh z/\lambda$, so $Q=CG$, $\tau_w=TG$, $\Delta p=LG$, and $\mathcal P=LGQ$. The fixed-gradient baseline is $G=1.0\ \mathrm{N/m^3}$. The optional fixed-flow baseline $Q=10^{-7}\ \mathrm{m^3/s}$ gives $G=0.688941996\ \mathrm{N/m^3}$, $\Delta p=0.068894200\ \mathrm{Pa}$, $\tau_w=8.894200\times10^{-4}\ \mathrm{Pa}$, and $\mathcal P=6.889420\times10^{-9}\ \mathrm W$.

Let $D_m=m\,\partial/\partial m$ with other inputs held fixed. The plate term contributes $dK_p=K_p(d\mu/\mu-2\,db/b)$. Differentiate $z=(H/2)\sqrt{K/\mu}$ before differentiating $C$ and $T$.

| Varied input $m$ | $D_m\log C$ | $D_m\log T$ |
|---|---:|---:|
| $\mu$ | −0.814823580 | 0.099999420 |
| $b$ | 2.481411359 | 0.799995364 |
| $\gamma$ | 0.277764630 | 0.149999131 |
| $\alpha_{\max}$ | −0.185176420 | −0.099999420 |
| $q_\alpha$ | 0.092588210 | 0.049999710 |
| $H$ | 1.148235801 | $5.795280\times10^{-6}$ |
| $L$ or $\rho$ | 0 | 0 |

At fixed $G$, $D_m\log Q=D_m\log C$ and $D_m\log\tau_w=D_m\log T$. At fixed $Q$, $D_m\log G=-D_m\log C$ and $D_m\log\tau_w=D_m\log T-D_m\log C$. Differentiate $\Delta p=LG$ and $\mathcal P=LGQ$ for drop and power. If $G$ itself varies, $Q$, $\tau_w$, and $\Delta p$ scale with $G$, and power scales with $G^2$.

Use one input at a time with $m_\pm=m_0(1\pm\epsilon)$ and $\epsilon=10^{-2},10^{-3},10^{-4},10^{-5}$. Compare $[J(m_+)-J(m_-)]/(2\epsilon)$ to $D_mJ$. Check the corrected Taylor remainder $|J(m_+)-J(m_0)-\epsilon D_mJ|$ for second-order scaling where numerical error is small. Compare tangent, adjoint, and differences on the same mesh before comparing with the continuum references under refinement.

For residual $\mathcal R(w,m)=0$, the scaled tangent solves $\mathcal R_wD_m w=-m_0\mathcal R_m$. The adjoint solves $\mathcal R_w^\mathsf T\lambda=J_w^\mathsf T$, giving $D_mJ=m_0(J_m-\lambda^\mathsf T\mathcal R_m)$. At fixed volumetric flow, include $G$ as a state unknown and impose the target-flow residual with the perturbed depth.

## Solver comparison log

| Rung | Record later | Solver measurements | Pass/fail |
|---|---|---|---|
| 1 | Velocity profile, section flow, signed wall shear, input power | | |
| 2 | Actual mesh sizes and triangle counts, three errors and observed orders | | |
| 3 | Signed flows, all force and power terms, balance residuals | | |
| 4 | Pulled-back field discrepancies and corresponding error norms | | |
| 5 | Jacobian-action error, residuals, exact-field errors, accepted step lengths | | |
| 6 | Control, parameter, response, tangent, adjoint, finite differences, Taylor remainders | | |

Solver tolerances and quantitative acceptance thresholds are pending.
