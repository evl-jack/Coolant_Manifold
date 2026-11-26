"""Independent mechanical energy quadratures for original and feasible states."""
import json
import sys
from pathlib import Path

import numpy as np

from .fem import quadratic_basis
from .project43_baseline import Geometry
from .project43_brinkman import create_problem


def audit(fem, ports, state, drag, scale):
    re = fem.inertia
    u, grad = fem.velocity_fields(state)
    sym = (grad + grad.swapaxes(-1, -2)) / 2
    speed2 = (u*u).sum(axis=-1)
    pressure = state[2*fem.n_u:][fem.triangles] @ fem.p_basis.T
    def integrate(field):
        return float(scale * np.sum(fem.weights * field))
    visc = integrate(2*(sym*sym).sum(axis=(-2,-1)))
    friction = integrate(drag*speed2)
    defect = integrate((pressure + re*speed2/2)*np.trace(grad, axis1=-2, axis2=-1))
    incident = {}
    for e, tri in enumerate(fem.triangles):
        for i,j in ((0,1),(1,2),(2,0)):
            incident[tuple(sorted((tri[i],tri[j])))] = e
    q,w = np.polynomial.legendre.leggauss(5)
    s,w = (q+1)/2,w/2
    terms = []
    for edges,n in zip(ports, ([-1,0],[1,0],[0,1])):
        term = np.zeros(3)
        for a,b,_ in edges:
            e = incident[tuple(sorted((a,b)))]
            tri = fem.triangles[e]
            bary = np.zeros((len(q),3))
            bary[:,np.flatnonzero(tri==a)[0]] = 1-s
            bary[:,np.flatnonzero(tri==b)[0]] = s
            basis, refgrad = quadratic_basis(bary)
            pts = fem.points[tri]
            g = refgrad @ np.linalg.inv((pts[1:]-pts[0]).T)
            vnodes = fem.velocity_dofs[e]
            nodal = np.column_stack((state[vnodes],state[vnodes+fem.n_u]))
            vel = basis @ nodal
            dv = np.einsum('qaj,ai->qij',g,nodal)
            tr = (dv+dv.swapaxes(-1,-2)) @ n
            vn = vel @ n
            pres = bary @ state[2*fem.n_u+tri]
            fields = np.stack((-pres*vn,-.5*re*(vel*vel).sum(axis=1)*vn,(vel*tr).sum(axis=1)))
            term += scale*np.linalg.norm(fem.points[a]-fem.points[b])*(fields @ w)
        terms.append(term.tolist())
    p,k,t = np.sum(terms,axis=0)
    return dict(units='W/m depth', ports_inlet_main_side_pressure_kinetic_viscous=terms,
                port_power=float(p+k), viscous_work=float(t), viscous_dissipation=visc,
                brinkman_dissipation=friction, weighted_divergence=defect,
                raw_closure=float(p+k+t-visc-friction),
                corrected_closure=float(p+k+t+defect-visc-friction))


def run(source, output):
    geom = Geometry()
    fem, base, fixed, values, ports, flux, _ = create_problem(geom,.5)
    mu,rho,height = .001002,998.2,.004
    scale = mu*(50*mu/(rho*height))**2
    alpha = 1e6*height**2/mu
    with np.load(source/'best.npz') as data:
        np.testing.assert_array_equal(data['triangles'],fem.triangles)
        np.testing.assert_allclose(data['vertices_mm'],fem.points*geom.height_mm)
        feasible = data['state'].copy()
        gamma = data['gamma'][:,None]
    result = {}
    for name,fraction in [('original',base),('feasible',gamma)]:
        drag = alpha*(1-fraction)**2
        if name=='original':
            state = None
            for fem.inertia in (0.,12.5,25.,37.5,50.):
                state = fem.solve(fixed,values,fem.load(),resistance=drag,initial=state).values
        else:
            state = feasible
        result[name] = audit(fem,ports,state,drag,scale)
        q = flux @ state
        result[name]['side_share'] = float(q[2]/q[0])
    with output.open('x') as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
    print(json.dumps(result))


if __name__=='__main__':
    run(Path(sys.argv[1]),Path(sys.argv[2]))
