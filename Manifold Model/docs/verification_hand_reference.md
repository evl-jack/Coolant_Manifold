# Six-rung hand reference — symbolic, before numerical verification

No numerical problem values are assigned here. The exact answers below are
references to derive on paper, not completed solver results. Convergence rates,
Newton histories, rotation errors, and gradient discrepancies require later
computed results; those cannot honestly be filled in by hand in advance.

Use distinct symbols: rho is physical density, gamma is the Brinkman design
field, z is the discrete state, and m is an independent input parameter.

## 1. Exact solutions: integrate a straight-channel equation

First distinguish two geometries. The old gap-resolved channel meshes the
distance between plates directly. The new top-down model averages over a
separate depth b. Do not use the old drag-free formula to verify a model with
plate drag active.

For a straight top-down strip, let x run along length L, y across in-plane
width H, and b be the unmeshed depth. Assume fully developed laminar flow
U=(u(y),0), constant viscosity, stationary side walls, and G=-dp/dx.
Use periodic ends plus body force G for a particularly clean reference.
Then convection is zero and the model reduces to

$$
\mu u''-cu+G=0,\qquad u(0)=u(H)=0,
\qquad c=\frac{12\mu}{b^2}+\alpha_0.
$$

alpha_0 is a prescribed spatially uniform Brinkman coefficient. For the
initial hollow interior alpha_0=0, but plate drag is still present.
Set kappa=sqrt(c/mu), solve the constant-coefficient ODE, and apply symmetry:

$$
u(y)=\frac{G}{c}\left[
1-\frac{\cosh(\kappa(y-H/2))}{\cosh(\kappa H/2)}\right].
$$

Integrating across the strip gives

$$
Q=b\int_0^H u\,dy
=\frac{bG}{c}\left[H-\frac{2}{\kappa}
\tanh\left(\frac{\kappa H}{2}\right)\right].
$$

The maximum is u(H/2), and average velocity is Q/(bH). Differentiate the
profile to obtain side-wall shear:

$$
|\tau_w|=\frac{\mu|G|\kappa}{c}
\tanh\left(\frac{\kappa H}{2}\right).
$$

For the separate drag-free case c=0, integrate mu*u''=-G twice:

$$
u(y)=\frac{G}{2\mu}y(H-y),\qquad
Q=\frac{bGH^3}{12\mu},\qquad u_{\max}=\frac{GH^2}{8\mu}.
$$

Use that expression directly at c=0 rather than substituting c=0 into the
hyperbolic formula. It is also the small-c limit.

For future 3D hollow circular-pipe verification, radial momentum gives

$$
\mu\frac1r\frac{d}{dr}\left(r\frac{du}{dr}\right)=-G,
\qquad u(R)=0,\quad u'(0)=0,
$$

$$
u(r)=\frac{G}{4\mu}(R^2-r^2),\qquad
Q=2\pi\int_0^R u(r)r\,dr=\frac{\pi R^4G}{8\mu}.
$$

This pipe reference assumes no design drag and resolves the wall in 3D;
it has no extra depth-averaged plate-drag term.

## 2. Convergence: prescribe a smooth field and derive its forcing

A polynomial channel profile can be represented exactly by P2 velocity and
is therefore unsuitable for measuring its general convergence rate.
Use a manufactured closed-box solution instead. With X=x/S, Y=y/S, and
A(t)=t^2(1-t)^2, prescribe

$$
u=U_0 A(X)A'(Y),\quad v=-U_0 A'(X)A(Y),\quad
p=P_0(X^3+Y^3-1/2).
$$

By hand, differentiate to show u_x+v_y=0, check u=v=0 at all walls, and
integrate p to show its mean is zero. The useful derivatives are

$$
A'=2t-6t^2+4t^3,\quad A''=2-12t+12t^2,\quad A'''=-12+24t.
$$

For linear Stokes plus constant drag, construct the force independently:

$$
\mathbf f=\nabla p-\mu\nabla^2\mathbf u+c\mathbf u.
$$

For example,

$$
f_x=\frac{3P_0X^2}{S}
-\frac{\mu U_0}{S^2}[A''(X)A'(Y)+A(X)A'''(Y)]
+cU_0A(X)A'(Y),
$$

$$
f_y=\frac{3P_0Y^2}{S}
+\frac{\mu U_0}{S^2}[A'''(X)A(Y)+A'(X)A''(Y)]
-cU_0A'(X)A(Y).
$$

Independent hand integrals for normalization are

$$
\|\mathbf u\|_{L^2}^2=\frac{2U_0^2S^2}{33075},\qquad
\|p\|_{L^2}^2=\frac{9P_0^2S^2}{56}.
$$

For smooth solutions on suitable meshes, expect P2/P1 velocity L2 error
O(h^3), pressure L2 error O(h^2), and velocity-gradient L2 error O(h^2).
After solving on successive meshes, calculate

$$
r=\frac{\log(E(h)/E(h/2))}{\log2}.
$$

On paper prepare the exact fields, forcing, norms, and expected rates.
Leave actual error and rate columns empty until numerical solves exist.

## 3. Global balance: integrate conservation laws

For the manifold, with outward unit normal n,

$$
Q_{in}=-b\int_{\Gamma_{in}}\mathbf U\cdot\mathbf n\,ds,
\qquad Q_i=b\int_{\Gamma_i}\mathbf U\cdot\mathbf n\,ds,
\qquad Q_{in}=\sum_i Q_i.
$$

In 3D replace b*ds with dA. Equal outlet flows are a performance target,
not a consequence of conservation or equal outlet pressures.

For the straight-strip reference, retain signed forces:

$$
F_{drive}=GLHb,
\quad F_{wall}=\frac{2bL\mu G\kappa}{c}\tanh(\kappa H/2),
\quad F_{drag}=cLQ,
\quad F_{drive}=F_{wall}+F_{drag}.
$$

Multiply its ODE by u and integrate by parts. Since u=0 at the walls,

$$
GLQ=bL\int_0^H[\mu(u')^2+cu^2],dy.
$$

This is the hand energy reference including plate and Brinkman losses.
The c=0 wall force follows directly from the parabolic profile.

For a general open 3D Navier-Stokes domain, do not equate pressure work to
viscous loss without the other terms. With no body force,

$$
-\int_{\partial\Omega}p\,\mathbf u\cdot\mathbf n\,dA
+\int_{\partial\Omega}\mathbf u\cdot(2\mu\varepsilon\mathbf n)\,dA
=\int_\Omega[2\mu\varepsilon:\varepsilon+\alpha|\mathbf u|^2],dV
+\frac{\rho}{2}\int_{\partial\Omega}|\mathbf u|^2\mathbf u\cdot\mathbf n\,dA.
$$

Here epsilon=sym(grad u). Boundary viscous work and net kinetic-energy
transport need not vanish. RANS has additional modeled energy transfers;
the laminar identity is not the complete RANS energy budget.

## 4. Invariance: transform the entire problem on paper

For a rotation matrix T with T^T*T=I and det(T)=1, define

$$
\mathbf x'=T\mathbf x,\quad
\mathbf u'(\mathbf x')=T\mathbf u(T^T\mathbf x'),\quad
p'(\mathbf x')=p(T^T\mathbf x').
$$

Rotate forces, boundary velocities, normals, and mesh coordinates too.
Scalars such as viscosity and alpha do not rotate. Derive

$$
\nabla'\mathbf u'=T(\nabla\mathbf u)T^T,
\qquad \mathbf u'\cdot\mathbf n'=\mathbf u\cdot\mathbf n.
$$

Therefore flow, pressure differences, and scalar power are unchanged.
For a 2D quarter turn, (x,y)->(-y,x) and (u,v)->(-v,u).
For the depth-averaged model, rotate in the modeled plane; a tilt of the
averaging direction is not the same physical 2D problem. Numerical errors
between transformed solutions are measured later.

## 5. Nonlinear benchmark: derive inertia and the Newton tangent

Use the manufactured fields from Rung 2 and add the exact inertial force:

$$
\mathbf f=\rho(\mathbf u\cdot\nabla)\mathbf u
+\nabla p-\mu\nabla^2\mathbf u+c\mathbf u.
$$

The acceleration components, derived independently on paper, are

$$
a_x=\frac{U_0^2}{S}A(X)A'(X)
\{[A'(Y)]^2-A(Y)A''(Y)\},
$$

$$
a_y=\frac{U_0^2}{S}A(Y)A'(Y)
\{[A'(X)]^2-A(X)A''(X)\}.
$$

The exact derivative of convection in a perturbation delta_u is

$$
DN[\mathbf u]\,\delta\mathbf u
=\rho[(\delta\mathbf u\cdot\nabla)\mathbf u
+(\mathbf u\cdot\nabla)\delta\mathbf u].
$$

Both terms are required. Later compare this Jacobian action with centered
differences of the full residual at a nonconverged state, and examine local
full-step Newton convergence above the roundoff floor. The expected local
behavior is ||e_{k+1}||=O(||e_k||^2), under the usual smoothness and
nonsingular-Jacobian conditions. Deriving it is not evidence that a particular
implementation achieves it. A RANS benchmark is a separate addition.

## 6. Gradients: start with a hand response derivative

For the separate drag-free straight channel at fixed G,

$$
Q=\frac{bGH^3}{12\mu},\quad
Q_G=Q/G,\quad Q_b=Q/b,\quad Q_H=3Q/H,\quad Q_\mu=-Q/\mu.
$$

These quotient expressions assume nonzero quantities; differentiate the
polynomial expression directly at G=0. If G is an independent input, Q_L=0.
If delta_p is held fixed instead, G=delta_p/L and Q_L=-Q/L. State what is held
fixed before differentiating.

For fixed Q in the same drag-free reference,

$$
P=\Delta p Q=\frac{12\mu LQ^2}{bH^3},\quad
P_\mu=P/\mu,\quad P_L=P/L,\quad P_Q=2P/Q,
\quad P_b=-P/b,\quad P_H=-3P/H.
$$

For the actual depth-averaged model, differentiate the c>0 exact expression
instead and include c_b=-24*mu/b^3 and c_mu=12/b^2 when alpha is independent.

For a general discrete response J(z(m),m), derive

$$
R(z,m)=0,\quad R_z z_m=-R_m,
\qquad R_z^T\lambda=J_z^T,
\qquad \frac{dJ}{dm}=J_m-\lambda^TR_m.
$$

The parameter derivative includes changed inlet values, outlet loads, mesh
coordinates, normal vectors, and integration measures when relevant.
Then perform, after implementation,

$$
D_{FD}=\frac{J(m+\epsilon d)-J(m-\epsilon d)}{2\epsilon},
\qquad D_{adj}=\nabla_mJ\cdot d,
$$

$$
E(\epsilon)=|J(m+\epsilon d)-J(m)-\epsilon\nabla_mJ\cdot d|
=O(\epsilon^2).
$$

Halving epsilon should reduce the corrected remainder by approximately four
until solver error or roundoff dominates. Use feasible directions at bounds.
Every perturbed evaluation must solve the changed PDE; holding velocity
fixed checks only a partial derivative. For field gradients specify whether
the returned gradient is a vector of coefficient derivatives or an L2
gradient (related by the mass matrix).

The Taylor-test procedure is documented in the
[dolfin-adjoint verification guide](https://dolfin-adjoint.github.io/dolfin-adjoint/documentation/verification.html).

## Report template

For each applicable scalar: test name, mesh, computed quantity, analytical
reference, absolute error, relative error, units, tolerance, and pass/fail.
Relative error is undefined at zero reference; use a stated dimensional
absolute tolerance instead. Convergence and invariance tests have expected
rates or relationships rather than a single exact scalar number.
