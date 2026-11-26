# Simultaneous baffle and full-topology diagnostics

Both models use the original envelope, ports and sharp T initialization,
nonlinear Navier–Stokes at Re = 50, an 83 mm² area cap and a 49.9–50.1%
side-flow band. Model A minimizes \(P/P_0\). Model B minimizes
\(P/P_0+0.1s\), where \(s=h/H\), \(H=4\) mm and \(0\le s\le0.9375\).
Both models use the same original sharp-initialized zero-height immersed
power \(P_0\). The hybrid zero-height evaluation records the change due to
the reserved strip before the bounded positive-height sweep and common
feasibility restoration.

## Filter and projection

The strip \(13\le x\le14\), \(0\le y\le4\) mm has pinned raw and
projected densities equal to one. The coarse fixture has 224 free topology
densities and one height coordinate. Let \(E\) inject free raw coordinates,
\(r_f\) be one in pinned cells and zero elsewhere, and \(D=EE^T\).
The full-domain Neumann filter retains the pinned cells in its stencil.
With cell mass \(M\), existing Helmholtz matrix \(A_f\) and projection
\(p_\beta\), the forward and transpose operations are

\[
 z=A_f^{-1}M(Ed+r_f),\quad
 \gamma=Dp_\beta(z)+r_f,\quad
 J_\gamma^T v=E^TM^TA_f^{-T}\operatorname{diag}(p_\beta'(z))Dv .
\]

The projected override zeros the pinned rows before the transpose solve.
The raw gradient is restricted to free coordinates after that solve.

## Residual and area

The clipped quintic equals \(S(t)=10t^3-15t^4+6t^5\) for \(0<t<1\)
and is constant zero and one outside. The baffle mask is
\(B=I_b(x)S((Hs-y)/\ell_s)\) with \(\ell_s=0.125\) mm and
\(B_s=I_b(x)(H/\ell_s)S'((Hs-y)/\ell_s)\).
Here \(S'=30t^2(1-t)^2\) inside and zero outside. The effective drag
coefficient is \((1-B)\alpha(\gamma)+B\alpha_{\max}\), where
\(\alpha=\alpha_{\max}q(1-\gamma)/(q+\gamma)\).
At fixed state \(u\), its contributions to the residual derivatives are

\[
 R_s[v]=\int B_s(\alpha_{\max}-\alpha(\gamma))u\cdot v,\quad
 R_d[\delta d;v]=\int(1-B)\alpha'(\gamma)
 (J_\gamma\delta d)u\cdot v,\quad
 \alpha'=-\frac{\alpha_{\max}q(q+1)}{(q+\gamma)^2}.
\]

The implementation multiplies by the existing nondimensional drag scale.
Within the strip, \(\alpha=0\) and \(J_\gamma\delta d=0\). The height
term uses the integrated baffle mask and the topology term uses the stock
element quadrature. Horizontal subcells split at element vertices,
\(y=h-\ell_s\) and \(y=h\). Three Gauss points across and six along height
integrate the P2 basis products times the quintic exactly, including the
linear section-width factor. Adjacent mask traces agree at the moving
limits, so internal boundary terms cancel.

The immersed fluid area integrates \((1-B)\gamma\) and includes the fixed
41 mm² of fluid. Its derivatives are \(A_s=-\int B_s\) and
\(A_d[\delta d]=\int J_\gamma\delta d\). The strip has width \(w=1\) mm.
Its smooth removed area is \(w\ell_s(t^6-3t^5+2.5t^4)\) at
\(t=h/\ell_s\le1\), and \(w(h-\ell_s/2)\) above that range.
At zero height the mask and its first derivative vanish. The independent
hard area and connectivity checks restore the exact reserved rectangle
to the \(\gamma=0.5\) reconstructed polygon before subtracting the exact
baffle. This preserves the strip even at its marching-squares corners.

## Hydraulic sensitivities

For \(F=P\) or \(F=f=Q_{\rm side}/Q_{\rm in}\), solve
\(R_U^T\lambda_F=F_U^T\) on the free state degrees of freedom.
The ports are fixed, so the total sensitivities are
\(-\lambda_F^TR_s\) and \(-\lambda_F^TR_d\).
The flow state gradient is
\((Q_{{\rm side},U}-f Q_{{\rm in},U})/Q_{\rm in}\).
Normalization divides power derivatives by the common \(P_0\).
Model B adds \(0.1\) to the derivative with respect to \(s\).

## Baffle-aware connectivity

For each graph segment, compact quintic gates \(a,b,c\) indicate
clearance lateral to the baffle, below the base, and above its tip.
The tip clearance uses the lowest point on the entire edge segment.
Each gate has a positive geometric margin.
Their union is \(m_e=1-(1-a)(1-b)(1-c)\) and
\(m_{e,h}=(1-a)(1-b)c_h\), where the tip clearance decreases at
unit rate with height. Every edge intersecting the hard rectangle
has zero clearance gate.

The new edge conductance is \(k_e=\epsilon_e+\kappa_e b_i b_j m_e\),
where \(b_i,b_j\) are the original projected-density gates. For the
minimized Dirichlet energy \(C_t\) to each outlet,

\[
 C_{t,s}=H\sum_e\kappa_e b_i b_j m_{e,h}(\phi_i-\phi_j)^2 .
\]

The density derivative multiplies each original graph-edge term by
\(m_e\) and then uses \(J_\gamma^T\). Free-potential derivatives cancel
by stationarity. A cut with no eligible edge has conductance bounded by
the original floor budget, strictly below the certificate threshold.
Both graph inequalities and the exact reconstructed polygons are checked.
Physical acceptance requires independent fitted-wall FEM, the existing
flow-split uncertainty checks and the original area and state tolerances.
The corrected Brinkman binary fixture remains a separate verification
task with its own thresholds.
