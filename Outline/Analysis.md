

# Anaylsis of Naviar Stokes w/ Brinkman penalty


## Overveiw

To summerize the problem We are looking at a straigtforward manifold with one inlet and 3 outlets. The model that we will be building is a 2D top down of a constant header depth. We will start with with steady, incompressiable, laminar Navier-Stokes flow. We will then use a depth averaged model to account for drag on the top and bottom plates.

- Inlet: Prescribe a velocity profile to initilize total flow $ Q_in $
- Walls : Apply no slip along the side walls
- Outlets : Give all three branches the same downstream resistance to start. The solver Determines $Q_1$ , $ Q_2 $ , $ Q_3 $. Equal resistance gives the branches the same downstream conditions, but their flows can still differ because pressure varies inside of the header.


## Analysis Model

### What am I modeling

#### Coolant Distribution System

I am modeling a coolant distribution manifold with n inlets and k outlets, I want to optimize geometry by placing a Brinkman Penalty on the governing Naviar-Stokes equation to enforce a drag force, which will model the walls inside of the manifold. By optimizing geomentry the goal is to minimize pumping power under the guide that the mass Flow of the Inlets are within a tolerance of eachother and sum to the total mass flow in.
- Using the Naviar-Stokes Equations to predict velocity and desnisty 
- Enforcing a Brinkman Penalty to model geometry and build density feilds
- Evaluating Preformance of optimized models under turbulant opperating conditions using RANS


### Governing Equations

#### Navier Stokes : Laminar Flow

$$
\rho(\mathbf{U}\cdot\nabla)\mathbf{U}
=
-\nabla p
+\mu\nabla^2\mathbf{U}
-\underbrace{\frac{12\mu}{b^2}\mathbf{U}}_{\text{top and bottom plate drag}}.
$$

$$
\nabla\cdot\mathbf{U}=0.
$$



#### Navier Stokes : Laminar Flow w/ Brinkman Penalization

$$
\rho(\mathbf{U}\cdot\nabla)\mathbf{U}
=
-\nabla p
+\mu\nabla^2\mathbf{U}
-\frac{12\mu}{b^2}\mathbf{U}
-\underbrace{\alpha(\gamma)\mathbf{U}}_{\text{Brinkman drag}}.
$$

$$
\nabla\cdot\mathbf{U}=0.
$$

The proposed interpolation for the Brinkman drag coefficient is

$$
\alpha(\gamma)
=
\alpha_{\max}\frac{1-\gamma}{1+q_\alpha\gamma},
\qquad
0\leq\gamma\leq1,
\qquad
q_\alpha>0.
$$

- For fixed geometry, any case where we would only use the Laminar Flow eqaution, we can set the Brinkman drag coeffient equal to zero and explicitly define the walls.

- will use the full equation by setting $ \gamma = 0 $ for any element restricting flow and $ \gamma = 1 $ for elements granting open flow.
- Must confirm that high drag suppresses velocity

#### Turbulant Flow unsing RANS

 mean momentum and continuity equations are

$$
\begin{aligned}
\rho(\overline{\mathbf{u}}\cdot\nabla)\overline{\mathbf{u}}
&=\nabla\cdot\boldsymbol{\sigma}_{\mathrm{RANS}},\\
\nabla\cdot\overline{\mathbf{u}}&=0.
\end{aligned}
$$

 mean stress and strain-rate tensors are

$$
\begin{aligned}
\boldsymbol{\sigma}_{\mathrm{RANS}}
&=-\left(\bar{p}+\frac{2}{3}\rho k\right)\mathbf{I}
+2(\mu+\mu_t)\mathbf{S},\\
\mathbf{S}
&=\frac{1}{2}\left[
\nabla\overline{\mathbf{u}}
+\left(\nabla\overline{\mathbf{u}}\right)^{\mathsf T}
\right].
\end{aligned}
$$

#### Menter k–omega SST Closure

Use the SST-2003 variant to calculate $k$, the specific dissipation rate
$\omega$, and the turbulent viscosity $\mu_t$.

$$
\rho\,\overline{\mathbf{u}}\cdot\nabla k
=
\widetilde{P}_k-\beta^*\rho k\omega
+\nabla\cdot\left[(\mu+\sigma_k\mu_t)\nabla k\right].
$$

$$
\begin{aligned}
\rho\,\overline{\mathbf{u}}\cdot\nabla\omega
&=
\frac{\gamma_\omega}{\nu_t}\widetilde{P}_k
-\beta\rho\omega^2
+\nabla\cdot\left[(\mu+\sigma_\omega\mu_t)\nabla\omega\right]\\
&\quad+
2(1-F_1)\frac{\rho\sigma_{\omega 2}}{\omega}
\,\nabla k\cdot\nabla\omega.
\end{aligned}
$$

The strain-rate magnitude, turbulent viscosity, and limited turbulence
production are

$$
\begin{aligned}
S&=\sqrt{2\mathbf{S}:\mathbf{S}},\\
\mu_t&=\frac{\rho a_1 k}{\max(a_1\omega,SF_2)},\\
\nu_t&=\frac{\mu_t}{\rho},\\
\widetilde{P}_k&=\min\left(\mu_t S^2,\;10\beta^*\rho k\omega\right).
\end{aligned}
$$

The blending functions $F_1$ and $F_2$ depend on the distance to the
nearest wall and on the turbulence variables. Use the published SST-2003
definitions for these functions and the model coefficients.
The SST coefficient $\gamma_\omega$ is separate from the Brinkman
design field $\gamma$.

[Menter SST model equations and coefficients](https://tmbwg.github.io/turbmodels/sst.html)

### Boundary Conditions and Flow Rates

Let $\mathbf{n}$ be the outward unit normal. In the 2D model, each port
is a line segment $\Gamma$ and $b$ is the manifold depth.

#### Inlet

Prescribe the inlet velocity profile $\mathbf{U}_{\mathrm{in}}$ to
deliver the chosen total volumetric flow rate.

$$
\mathbf{U}=\mathbf{U}_{\mathrm{in}}
\quad\text{on }\Gamma_{\mathrm{in}},
\qquad
Q_{\mathrm{in}}
=-b\int_{\Gamma_{\mathrm{in}}}
\mathbf{U}_{\mathrm{in}}\cdot\mathbf{n}\,ds.
$$

The minus sign makes $Q_{\mathrm{in}}$ positive because incoming flow
points opposite the outward normal. For RANS, also specify inlet
values for $k$ and $\omega$.

#### Stationary Walls

Apply no slip along explicit stationary walls.

$$
\mathbf{U}=\mathbf{0}
\quad\text{on }\Gamma_w
\qquad\text{(2D laminar model)}.
$$

$$
\overline{\mathbf{u}}=\mathbf{0},
\qquad
k=0
\quad\text{on }\Gamma_w
\qquad\text{(wall-resolved RANS)}.
$$

Use the SST near-wall condition for $\omega$ in the RANS model.

#### Outlets

The volume flow through outlet $i$ is

$$
Q_i=b\int_{\Gamma_i}\mathbf{U}\cdot\mathbf{n}\,ds.
$$

In the depth-resolved RANS model, integrate over the actual port areas
instead.

$$
Q_{\mathrm{in}}
=-\int_{A_{\mathrm{in}}}\overline{\mathbf{u}}\cdot\mathbf{n}\,dA,
\qquad
Q_i=\int_{A_i}\overline{\mathbf{u}}\cdot\mathbf{n}\,dA.
$$

Represent the cooling circuit beyond outlet $i$ by a pressure–flow
relationship $\Delta p_{\mathrm{down},i}(Q_i)$.

$$
p_{\mathrm{load},i}
=p_{\mathrm{ref}}+\Delta p_{\mathrm{down},i}(Q_i),
\qquad
\boldsymbol{\sigma}\mathbf{n}
=-p_{\mathrm{load},i}\mathbf{n}
\quad\text{on outlet }i.
$$

For laminar flow, the stress used in this condition is

$$
\boldsymbol{\sigma}
=-p\mathbf{I}
+\mu\left[
\nabla\mathbf{U}+(\nabla\mathbf{U})^{\mathsf T}
\right].
$$

For turbulent flow, use $\boldsymbol{\sigma}_{\mathrm{RANS}}$.
For a straight round downstream pipe with fully developed laminar flow,
the pressure loss is

$$
\Delta p_{\mathrm{down},i}=R_iQ_i,
\qquad
R_i=\frac{8\mu L_i}{\pi r_i^4},
$$

where $L_i$ and $r_i$ are the pipe length and radius. For a turbulent
downstream circuit, use a pressure–flow relationship appropriate to
that operating condition.

### Reported Quantities

#### Mass Balence

For a manifold with $N$ outlets, the outlet flows should sum to the
inlet flow. The normalized mass-balance error is

$$
E_{\mathrm{mass}}
=\frac{\left|Q_{\mathrm{in}}-\sum_{i=1}^{N}Q_i\right|}
{Q_{\mathrm{in}}}.
$$

#### Outlet Flow Uniformity

For branches with equal cooling demand, the target flow is

$$
Q_i^\star=\frac{Q_{\mathrm{in}}}{N}.
$$

The largest fractional deviation from the branch targets is

$$
E_Q=\max_{i=1,\ldots,N}
\frac{|Q_i-Q_i^\star|}{Q_i^\star}.
$$

#### Pressure Drop and Hydraulic Pumping Power

Calculate the flow-weighted inlet pressure. In the 2D model it is

$$
p_{\mathrm{in,fw}}
=\frac{b}{Q_{\mathrm{in}}}
\int_{\Gamma_{\mathrm{in}}}
p\,(-\mathbf{U}\cdot\mathbf{n})\,ds.
$$

For the depth-resolved RANS model, use the inlet area and mean fields.

$$
p_{\mathrm{in,fw}}
=\frac{1}{Q_{\mathrm{in}}}
\int_{A_{\mathrm{in}}}
\bar{p}\,(-\overline{\mathbf{u}}\cdot\mathbf{n})\,dA.
$$

The pressure drop to the common downstream reference and the
corresponding pressure-based hydraulic power are

$$
\Delta p=p_{\mathrm{in,fw}}-p_{\mathrm{ref}},
\qquad
P_{\mathrm{hyd}}=Q_{\mathrm{in}}\Delta p.
$$



## Verifications

### Exact solutions

| Quantity | Symbol | Value |
|---|---|---:|
| Gap between plates | $H$ | $0.001\ \mathrm{m}$ |
| Channel length | $L$ | $0.10\ \mathrm{m}$ |
| Channel width into the page | $W$ | $0.10\ \mathrm{m}$ |
| Driving pressure drop | $\Delta p$ | $1.2\ \mathrm{Pa}$ |
| Dynamic viscosity | $\mu$ | $0.001\ \mathrm{Pa\,s}$ |

Side-wall effects are neglected. The coordinate $x$ runs downstream,
and $y$ runs from the lower plate at $y=0$ to the upper plate at $y=H$.

Because the velocity profile is unchanged along the channel, convective
acceleration is zero. The downstream momentum equation reduces to

$$
\mu\frac{d^2u}{dy^2}
=\frac{dp}{dx}
=-\frac{\Delta p}{L}
=-12\ \mathrm{Pa/m}.
$$

Integrating twice gives

$$
u(y)=-\frac{\Delta p}{2\mu L}y^2+C_1y+C_2.
$$

No slip gives $u(0)=u(H)=0$. Therefore $C_2=0$ and
$C_1=\Delta p H/(2\mu L)$, yielding the exact velocity profile

$$
u(y)=\frac{\Delta p}{2\mu L}y(H-y).
$$

The maximum velocity occurs halfway between the plates:

$$
u_{\max}=u(H/2)
=\frac{\Delta p H^2}{8\mu L}
=\frac{1.2(0.001)^2}{8(0.001)(0.10)}
=0.0015\ \mathrm{m/s}.
$$

Integrating the velocity across the gap gives the volume flow rate:

$$
\begin{aligned}
Q
&=W\int_0^H u(y)\,dy\\
&=\frac{W\Delta p}{2\mu L}
  \left[\frac{Hy^2}{2}-\frac{y^3}{3}\right]_0^H\\
&=\frac{W\Delta p H^3}{12\mu L}\\
&=\frac{(0.10)(1.2)(0.001)^3}{12(0.001)(0.10)}
=1.0\times10^{-7}\ \mathrm{m^3/s}
=6.0\ \mathrm{mL/min}.
\end{aligned}
$$

The average velocity is

$$
\bar u=\frac{Q}{WH}
=\frac{1.0\times10^{-7}}{(0.10)(0.001)}
=0.0010\ \mathrm{m/s}.
$$

The shear-stress magnitude at either wall is

$$
|\tau_w|
=\left|\mu\frac{du}{dy}\right|_{\mathrm{wall}}
=\frac{\Delta p H}{2L}
=\frac{1.2(0.001)}{2(0.10)}
=0.006\ \mathrm{Pa}.
$$

The forward pressure force and the combined backward drag from both plates
have equal magnitudes:

$$
F_{\mathrm{pressure}}=\Delta p\,WH
=1.2\times10^{-4}\ \mathrm{N},
\qquad
F_{\mathrm{drag}}=2|\tau_w|WL
=1.2\times10^{-4}\ \mathrm{N}.
$$

Impose no slip on the plates and periodic conditions at the two ends.
Drive the flow with a uniform downstream force per unit volume
$G=\Delta p/L=12\ \mathrm{N/m^3}$, equivalent to the specified pressure
gradient. Use this force as the driving input so the solver predicts
the velocity and flow rate independently.

| Quantity | Exact reference | Computed result | Relative error |
|---|---:|---:|---:|
| Centerline velocity | $1.5\ \mathrm{mm/s}$ | — | — |
| Average velocity | $1.0\ \mathrm{mm/s}$ | — | — |
| Volume flow rate | $6.0\ \mathrm{mL/min}$ | — | — |
| Wall shear-stress magnitude | $0.006\ \mathrm{Pa}$ | — | — |
| Total wall drag magnitude | $1.2\times10^{-4}\ \mathrm{N}$ | — | — |

For each reported quantity $a$, calculate

$$
\text{Relative error}=
\frac{|a_{\mathrm{computed}}-a_{\mathrm{exact}}|}
     {|a_{\mathrm{exact}}|}\times100\%.
$$


### Convergence Rate

The straight-channel velocity profile is quadratic, so quadratic velocity
elements may reproduce it almost exactly on a coarse mesh. To measure
convergence, use a manufactured solution with higher-order spatial
variation and compare the computed fields against the exact fields.

#### Geometry and Exact Solution

Use a square domain $\Omega=[0,S]\times[0,S]$, where $S$ is the side
length. Define dimensionless coordinates and a polynomial function:

$$
X=\frac{x}{S},
\qquad
Y=\frac{y}{S},
\qquad
A(z)=z^2(1-z)^2.
$$

The derivatives needed for the forcing are

$$
\begin{aligned}
A'(z)&=2z-6z^2+4z^3,\\
A''(z)&=2-12z+12z^2,\\
A'''(z)&=-12+24z.
\end{aligned}
$$

Choose the exact velocity and pressure fields as

$$
\begin{aligned}
u_{\mathrm{exact}}&=U_0A(X)A'(Y),\\
v_{\mathrm{exact}}&=-U_0A'(X)A(Y),\\
p_{\mathrm{exact}}&=P_0\left(X^3+Y^3-\frac12\right).
\end{aligned}
$$

Here $U_0$ and $P_0$ are fixed velocity and pressure scales. The velocity
satisfies incompressibility exactly:

$$
\frac{\partial u_{\mathrm{exact}}}{\partial x}
+\frac{\partial v_{\mathrm{exact}}}{\partial y}
=
\frac{U_0}{S}A'(X)A'(Y)
-\frac{U_0}{S}A'(X)A'(Y)
=0.
$$

Since $A(0)=A(1)=A'(0)=A'(1)=0$, both velocity components vanish
along all four walls. The exact pressure has zero mean over the square.

#### Governing Equation and Manufactured Forcing

First test the linear Stokes problem:

$$
-\nabla p+\mu\nabla^2\mathbf{u}-c\mathbf{u}+\mathbf{f}
=\mathbf{0},
\qquad
\nabla\cdot\mathbf{u}=0.
$$

Use $c=0$ for the basic Stokes check. Repeat the test with
$c=12\mu/b^2$ to verify the depth-averaged plate-drag term.

Calculate the forcing from the exact fields before implementing the
solver:

$$
\mathbf{f}
=\nabla p_{\mathrm{exact}}
-\mu\nabla^2\mathbf{u}_{\mathrm{exact}}
+c\mathbf{u}_{\mathrm{exact}}.
$$

The resulting force components are

$$
\begin{aligned}
f_x
&=\frac{3P_0X^2}{S}
-\frac{\mu U_0}{S^2}
\left[A''(X)A'(Y)+A(X)A'''(Y)\right]
+cU_0A(X)A'(Y),\\
f_y
&=\frac{3P_0Y^2}{S}
+\frac{\mu U_0}{S^2}
\left[A'''(X)A(Y)+A'(X)A''(Y)\right]
-cU_0A'(X)A(Y).
\end{aligned}
$$

Impose no slip on all four walls and a zero-mean pressure condition:

$$
\mathbf{u}=\mathbf{0}\quad\text{on }\partial\Omega,
\qquad
\int_\Omega p\,dA=0.
$$

#### Mesh Refinement and Error Measurement

Use triangular $P_2/P_1$ Taylor–Hood elements. Solve the same problem
on progressively finer meshes with characteristic sizes
$h,\ h/2,\ h/4,\ldots$. Keep the domain, forcing, material properties,
and boundary conditions fixed.

Let $\mathbf{u}_h=(u_h,v_h)$ and $p_h$ be the computed fields.
Measure the absolute velocity and pressure errors using the $L^2$ norm:

$$
E_u(h)
=
\sqrt{
\int_\Omega
\left[
(u_h-u_{\mathrm{exact}})^2
+(v_h-v_{\mathrm{exact}})^2
\right]\,dA
}.
$$

$$
E_p(h)
=
\sqrt{\int_\Omega(p_h-p_{\mathrm{exact}})^2\,dA}.
$$

Normalize by the corresponding exact-field norms to report relative
percentage errors:

$$
e_u(h)=100\frac{E_u(h)}{\|\mathbf{u}_{\mathrm{exact}}\|_{L^2(\Omega)}},
\qquad
e_p(h)=100\frac{E_p(h)}{\|p_{\mathrm{exact}}\|_{L^2(\Omega)}}.
$$

For either error measure, calculate the observed convergence rate:

$$
r=\frac{\log(E(h)/E(h/2))}{\log 2}.
$$

#### Expected Results

For this smooth solution on sufficiently fine meshes, the expected
$P_2/P_1$ convergence behavior is

| Error measure | Expected behavior | Expected rate |
|---|---|---:|
| Velocity $L^2$ error | $E_u\propto h^3$ | $r_u\approx3$ |
| Pressure $L^2$ error | $E_p\propto h^2$ | $r_p\approx2$ |

Record the measured results without assigning them in advance:

| Mesh size | Element count | Velocity error (%) | Velocity rate | Pressure error (%) | Pressure rate |
|---|---:|---:|---:|---:|---:|
| $h$ | — | — | — | — | — |
| $h/2$ | — | — | — | — | — |
| $h/4$ | — | — | — | — | — |
| $h/8$ | — | — | — | — | — |

Plot the errors against $h$ on logarithmic axes and compare the slopes
with the expected rates. Use sufficiently accurate numerical integration
and solve the linear system tightly enough that mesh error dominates.

This rung verifies the Stokes discretization first. The nonlinear
benchmark will add convective inertia to the equations and update the
manufactured forcing accordingly.

### Global Balance

Calculate the inlet and outlet flows directly from the computed velocity field and verify that total inflow equals total outflow. Report the normalized imbalance and repeat the check as the mesh is refined. Unequal branch flows are acceptable for this check as long as the total balances.

For the straight-channel benchmark, integrate the computed wall shear and compare the total wall drag with the known driving force. When plate drag or Brinkman drag is active, include those contributions in the force balance. Establish acceptance tolerances before evaluating the solver.

### Invariance

Rotate the entire manufactured problem and its mesh by 90 degrees, including the exact velocity field, applied force, and boundary conditions. Verify that velocity vectors rotate with the geometry while pressure, speed, and integral performance quantities agree at corresponding locations, using the same pressure reference. Repeat after mesh refinement and report normalized differences.

### Nonlinear Benchmark

Extend the manufactured square-flow test to include convective inertia. Derive the forcing from the prescribed exact fields with the inertial term included, then compare the computed velocity and pressure with those fields.

Check the Newton Jacobian against finite differences at a nonconverged state. Record the residual history and look for quadratic convergence near the discrete solution when full Newton steps are accepted. Repeat the mesh study to distinguish nonlinear-solver error from discretization error.

Before applying RANS to the manifold, separately verify the selected SST implementation using a documented benchmark and assess pressure loss and near-wall mesh sensitivity on a suitable turbulent-flow reference case. Record whether the comparison uses an exact solution, reference computation, or experimental data.

### Gradient Check

Begin with the straight-channel flow-rate sensitivity derived by hand. Once design-dependent responses are implemented, compare each adjoint directional derivative with central finite differences for every objective, constraint, and kind of design variable.

Repeat over a range of perturbation sizes and several directions. Perform a Taylor test by subtracting the first-order adjoint prediction from the actual response change and check for second-order remainder behavior before numerical noise dominates. Include derivatives through any filtering, drag interpolation, and geometry changes used by the model. Record discrepancies and solver tolerances.
