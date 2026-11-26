"""Minimize net advective mechanical energy entering at the physical ports."""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from .project43_baseline import Geometry
from .project43_brinkman import create_problem


def port_power(fem, ports, state, re):
    # P2 trace for velocity and P1 for pressure. Four Gauss points integrate
    # the degree-six kinetic flux exactly along each straight port edge.
    q, w = np.polynomial.legendre.leggauss(4)
    t, w = (q+1)/2, w/2
    b = np.stack(((1-t)*(1-2*t), t*(2*t-1), 4*t*(1-t)), axis=1)
    p = np.stack((1-t, t), axis=1)
    grad = np.zeros(fem.size)
    powers = []
    for edges, normal in zip(ports, ([1,0], [-1,0], [0,-1])):
        normal = np.asarray(normal)
        weight = np.linalg.norm(fem.points[edges[:,1]]-fem.points[edges[:,0]], axis=1)[:,None]*w
        uv = np.stack([state[edges+i*fem.n_u]@b.T for i in range(2)], axis=-1)
        pn = edges[:,:2]+2*fem.n_u
        head = state[pn]@p.T + .5*re*np.sum(uv**2, axis=-1)
        un = uv@normal
        powers.append(float(np.sum(weight*un*head)))
        np.add.at(grad, pn, (weight*un)@p)
        for i in range(2):
            np.add.at(grad, edges+i*fem.n_u, (weight*(normal[i]*head+un*re*uv[:,:,i]))@b)
    return powers, grad


def run(source, out):
    ref = json.loads((source/'report.json').read_text())['best']
    initial = np.array(ref['design'])
    out.mkdir()
    g = Geometry()
    fem, base, fixed, prescribed, ports, flux, _ = create_problem(g,.5)
    x, y = (fem.points[fem.triangles].mean(axis=1)*g.height_mm).T
    throttle = (x>13)&(x<14)&(y<4)
    shoulders = (abs(x-10)>.5)&(abs(x-10)<1.5)&(y>4)&(y<6)
    direction = np.stack((-throttle.astype(float), shoulders.astype(float)),axis=-1)[:,None,:]
    area_grad = np.einsum('eq,eld->d',fem.weights,direction)*g.height_mm**2
    np.testing.assert_allclose(area_grad,[-4,4])
    re, mu, rho, h = 50., .001002, 998.2, .004
    speed = re*mu/(rho*h)
    alpha = 1e6*h*h/mu
    scale = mu*speed**2
    fem.inertia = re
    load = fem.load()
    with np.load(source/'best.npz') as saved:
        warm = saved['state'].copy()
        np.testing.assert_array_equal(saved['triangles'],fem.triangles)
        np.testing.assert_allclose(saved['gamma'],(base+direction@initial)[:,0], atol=1e-12)
    history, cache = [], {}
    def evaluate(z):
        key = tuple(z)
        if key not in cache:
            gamma = base+direction@z
            drag = alpha*(1-gamma)**2
            state = fem.solve(fixed,prescribed,load,resistance=drag,initial=warm)
            qw = flux@state.values
            imbalance = abs(qw[0]-qw[1:].sum())/qw[0]
            assert np.min(qw)>0 and imbalance<1e-9
            share = qw[2]/qw[0]
            powers, grad = port_power(fem,ports,state.values,re)
            total = sum(powers)
            response = np.stack([grad,(flux[2]-share*flux[0])/qw[0]])
            jac = fem.drag_adjoint(state,response,-2*alpha*(1-gamma[:,:,None])*direction)
            area = 83+area_grad@z
            record = dict(design=list(z), fluid_area_mm2=area, side_share=share,
                          power_W_per_m=scale*total, signed_port_powers_W_per_m=(scale*np.array(powers)).tolist(),
                          relative_mass_imbalance=imbalance, residual=state.residual_history[-1])
            history.append(record)
            cache[key] = total,share,jac,record
            (out/'history.json').write_text(json.dumps(history,indent=2))
        return cache[key]
    p0,s0,_,baseline = evaluate(initial)
    np.testing.assert_allclose(s0, ref['side_share'], atol=1e-8,rtol=0)
    assert .499 <= s0 <= .501 and baseline['fluid_area_mm2'] <= 83+1e-9 and p0>0
    def constraints(z):
        share = evaluate(z)[1]
        return np.array([share-.499,.501-share,-area_grad@z])
    def constraint_jac(z):
        ds = evaluate(z)[2][1]
        return np.stack([ds,-ds,-area_grad])
    result = minimize(lambda z:(evaluate(z)[0]/p0,evaluate(z)[2][0]/p0),
                      initial, method='SLSQP',jac=True,bounds=[(0,1)]*2,
                      constraints=[dict(type='ineq', fun=constraints,jac=constraint_jac)],
                      options=dict(maxiter=50,ftol=1e-12))
    feasible = [r for r in history if .499<=r['side_share']<=.501 and r['fluid_area_mm2']<=83+1e-9]
    best = min(feasible,key=lambda r:r['power_W_per_m'])
    (out/'report.json').write_text(json.dumps(dict(source=str(source),starting=baseline,best=best,
        objective='Net signed port flux of static pressure plus kinetic energy, W/m of depth',
        budget_mm2=83.,fractional_power_reduction=1-best['power_W_per_m']/baseline['power_W_per_m'],
        optimizer_success=bool(result.success),termination=str(result.message)),indent=2))


if __name__=='__main__':
    run(Path(sys.argv[1]),Path(sys.argv[2]))
