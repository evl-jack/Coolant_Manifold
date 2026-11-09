Project 4 asks for a real engineering problem bound by a PDE.


A cold plate \ Flow manifold problem has been selected. 

**Problem**: A flow manifold is an instrument used to distribute fluid or gas from one inlet into multiple chanels. Commenly used in cooling systems to distribute a coolent between multiple cooling elements. The measurable benchmarks for a flow manifold are as follows

Flow Distribution Uninifomity: How evenly the fluid or gas is distributed between outlets. A proper manifold seeks to achive below 5% differnce

Pressure Drop : Pressure lost between the inlet and outlet. Due to friction, turbulance, and restrictive geometry.

Discharge Ratio ( $$C_d$$ ) : the Ratio of actual flow rate to ideal flow rate

Volumetric Flow Rate: the total volume of mass flowing through the manifold per unit time.

**Summery**

-- For this Model I will be observing these preformance metrics first before thinking of introducing structural or mechanical benchmarks.

For this model the goal would be to optimize these metrics. Under certine to be defined assumptions. In a direct to chip cooling system there is usally three main parts. The cold plate, some distribution\collection manifold or CDU(Coolant Distribution Unit), and the heat exchanger system. I will assume for the problem that initially I will optimize a **** server racks and systemically increase the call for coolant channels as a means of testing the optimizer.

In an Ideal coolant system We want to maintain laminar flow when possible and minimize turbulance. The addition of turbulance creates a number of setbacks to the system.

Turbulance can cause pressure drops in the manifold making optimal design of junctions and shaping important.

Pressure drops means will be supported by a pump meaning, more pressure drop more pumping power is required adding to the overall cost.

**models----------------------------------------------------------------------------------------**

There is two types of cooling manifolds for a standard data center. There is the CDU to server manifold and there is the in-rack manifolds. the CDU manifold delivers coolant to the server where another manifold distributes the coolant directly to the chip.

Coolant:

• PG25 (25% Propylene Glycol + Water): Endorsed by industry standards bodies including ASHRAE TC 9.9 and the Open Compute Project (OCP). It provides an ideal balance of freeze protection, anti-corrosion properties, microbial control, and safety.


• Ethylene Glycol (EG): Offers slightly better thermal conductivity and lower viscosity for high-density AI clusters, but is used less frequently due to its acute toxicity.

**Test 1**
First I will look as a minimal case to build a base model to hopefully be built appon. I will first observe the fluid dynamics of CDU manifold with one inlet and three outlets. Prescribe a fixed total inlet flow \(Q\) and give the three outlets equal downstream resistance. Record the pressure drop and the flow through each outlet. Then vary one or two header dimensions to reduce pumping power while keeping each outlet’s flow within 5% of its target, \(Q/3\).



**Naviar Stokes & Birkman Penalization and fluid dynamic**

How can we use Naviar Stokes to govern and model 

From Fluid Dynamics the Reynolds number

$$ R_e = ( \rho * V  * D ) /  \mu $$
 or 
$$ R_e = ( V  * D ) /  v $$

for 

Re < 2300 : Laminar Flow 

2300 < Re < 4000 : Transional Flow

Re > 4000 : Turbulant Flow

Naviar stokes : ---------------------------------------------------------------------------------------

*Momentum Equation*: Newtons second law applied to the fluid. Tracks of the fluids velocity and momentum are changed over ttime by balencing fluids accelerations against all internal and external forces acting on it.

$$
\rho(\mathbf{u}\cdot\nabla)\mathbf{u}
= -\nabla p
+ \mu\nabla^2\mathbf{u}
- \alpha\bigl(\gamma(x,y)\bigr)\mathbf{u}.
$$

At the channel walls the velocity is 0

*Momentum Equation with Brinkman penalization*: The brinkman penailization is a drag force applied to the right side of the equation opposing the fluids motion acting as solid wall. This gives us a mean of controlling the flow of the mesh

$$
\rho(\mathbf{u}\cdot\nabla)\mathbf{u}
= -\nabla p
+ \mu\nabla^2\mathbf{u}
- \alpha\bigl(\gamma(x,y)\bigr)\mathbf{u}.
$$

With regular naivar stokes the wall velocity is set to 0. With Brinkman penilzation the drag coeffient \alpha is made vary large driving the velocity to 0.

*Conservation of Mass*: The strict constraint imposing the rule that fluid can not be created or destroyed

$$
\nabla\cdot\mathbf{u}=0.
$$

the conservation of mass condition remains unchanged with the addition of the brinkman penilization

----------------------------------------------------------------
 
*Birkman Density Feild*

The Birkman Denisty field is a map that tells the flow solver where coolant can move freely and where it should be blocked. For this problem the optimizer changes the map to create channels inside the header.

The feild is just a value that veries from place to place. A value \gamma describes the design. 

\gamme = 1  : an open passage that coolant flows through

\gamma = 0 a nearly blocked region that represents solid material

0 < \gamma < 1 : a region decided by the optimizer to make changes  

-Have a manifold split in half horizontaly and look at it face down. we see a 2D area. We can define this area a as *f_0*. This area is made up of a mesh with N amount of elements. each element has a value $\gamma_e$ wich defines the amount of area that cell takes up. We also define a value *fmax* which is the constraint to our optimizer telling it what the maximum amount of space the cells can take up. This builds our feild for a single plane (x,y,z) where we hold z at a constant value. to take it further we could extrude the model to show up the feild at each plane without hold z constant.

*Turbulance --RANS(Reynolds averaged Navier-Stokes)*

- in turbulant flow velocity fluctuates over time. RANS splits the instantaneouse velocity into an average and a fluctuation

$$
\mathbf{u}(x,y,t)
=
\underbrace{\overline{\mathbf{u}}(x,y)}_{\text{mean velocity}}
+
\underbrace{\mathbf{u}'(x,y,t)}_{\text{turbulent fluctuation}}.
$$

- When averaging the Navier-Stokes momentum equation, the nonlinear velocity term leaves an extra stress from the fluctuations. This is called the **Reynaolds Stress**

$$
\rho(\overline{\mathbf{u}}\cdot\nabla)\overline{\mathbf{u}}
=
-\nabla\bar{p}
+\mu\nabla^2\overline{\mathbf{u}}
-\underbrace{
\nabla\cdot\left(\rho\,\overline{\mathbf{u}'\otimes\mathbf{u}'}\right)
}_{\text{effect of turbulent fluctuations}}.
$$

- We can use $ \kappa $ - $ \omega $ SST to solve two new equations 

$$
\begin{aligned}
\text{Transport of } k
&= \text{production} - \text{dissipation} + \text{diffusion},\\
\text{Transport of } \omega
&= \text{production} - \text{destruction} + \text{diffusion}
   + \text{cross-diffusion}.
\end{aligned}
$$

- $k$ is the turbulent kinetic energy per unit mass.
- $\omega$ is the specific dissipation rate.
- $\mu_t$ is the turbulent viscosity calculated from $k$, $\omega$,
  and the mean flow. It is used to model the Reynolds stress.

To see full equation: (https://tmbwg.github.io/turbmodels/sst.html)

*USES*

We will use RANS on our Navier-Stokes-Brinkman model to investigate how it preforms under turbulant operating conditions. Versus using RANS to create a turbulant optimized brinkman model. 
- Using RANS is difficult becasue with brinkman we see the walls as areas with high drag within the mesh. To implement RANS into this model we would have to find a way to treat turbulance in these conditions in a constant way, which is somting turbulant flow doesnt do.

 **Optimization**------------------------------------------------------------------------------------------------------------------------------------------------------

Before we build the optimizer first we need to see if our proposed model passes the given verification ladder!!!!!


**Design Variables**


Flow rate per branch:

Flow variation/ disterbence uniformity:

Pressure Drop:

Working Pressure:

Burst Pressure :

Manifold Cross section :

Header size and Diameter :

Termination style :




**Resources**

https://www.opencompute.org/documents/ocp-white-paper-rack-manifold-requirements-and-qualification-v3-pdf/


