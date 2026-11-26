"""Area-preserving reconstruction of two relaxed controls into impermeable walls."""
import json
from pathlib import Path
import sys
import numpy as np
from .fem import TaylorHood
from .junction import boundary_data, strip_triangles
from .fillet import pressure_weights
from .project43_power import port_power


def axis(cuts, h):
    cuts = np.unique(cuts)
    return np.r_[np.concatenate([np.linspace(a,b,max(1,int(np.ceil((b-a)/h)))+1)[:-1]
                                 for a,b in zip(cuts[:-1],cuts[1:])]),cuts[-1]]


def solve(xaxis,yaxis,throat,shoulder):
    x,y = np.meshgrid(xaxis,yaxis)
    points = np.column_stack((x.ravel(),y.ravel()))
    tris = strip_triangles(np.arange(x.size).reshape(x.shape))
    cx,cy = points[tris].mean(axis=1).T
    main = (cy<4) & ~((cx>13)&(cx<14)&(cy>throat))
    stem = (cy>4) & (abs(cx-10)<(.5+shoulder*(cy<6)))
    ids,inverse = np.unique(tris[main|stem],return_inverse=True)
    fem = TaylorHood(points[ids]/4,inverse.reshape(-1,3))
    fixed,prescribed,ports,flux = boundary_data(fem)
    widths = [4*np.linalg.norm(fem.points[e[:,1]]-fem.points[e[:,0]],axis=1).sum() for e in ports]
    np.testing.assert_allclose(widths,[4,4,1])
    area = fem.weights.sum()*16
    np.testing.assert_allclose(area,83-(4-throat)+4*shoulder,atol=1e-10)
    state = None
    for re in (0,25,50):
        fem.inertia = re
        state = fem.solve(fixed,prescribed,fem.load(),initial=state.values if state else None)
    inlet,main,side = flux@state.values
    balance = abs(inlet-main-side)/inlet
    assert balance<1e-9 and min(main,side)>0
    pin,pout = pressure_weights(fem,ports)
    mu,rho,h = .001002,998.2,.004
    speed = 50*mu/(rho*h)
    loss = mu*speed/h*(pin@state.values-state.values@(pout@state.values)/(main+side))
    powers,_ = port_power(fem,ports,state.values,50)
    return dict(triangles=len(fem.triangles), area_mm2=area, port_widths_mm=widths,
                side_share=side/inlet, relative_mass_imbalance=balance, residual=state.residual_history[-1],
                pressure_loss_Pa=loss, pumping_power_W_per_m=sum(powers)*mu*speed**2)


def run(source,out):
    ref = json.loads((source/'report.json').read_text())['best']
    throttle, shoulder = ref['design']
    assert 0<=throttle<1 and 0<=shoulder<=1
    throat = 4*(1-throttle)
    xa = axis([0,9.5-shoulder,9.5,10.5,10.5+shoulder,13,14,20],.25)
    ya = axis([0,throat/4,throat/2,3*throat/4,throat,4,6,7],.25)
    out.mkdir()
    original = solve(xa,ya,4.,0.)
    candidate = solve(xa,ya,throat,shoulder)
    np.testing.assert_allclose(candidate['area_mm2'],ref['fluid_area_mm2'],atol=1e-8)
    report = dict(source=str(source),relaxed=ref,throttle_open_height_mm=throat,
                  branch_widening_each_side_mm=shoulder,original=original,candidate=candidate,
                  mapping='Top-mounted rectangular baffle from x=13 to 14; widen both branch sides over y=4 to 6. Sharp corners, no-slip fitted walls.',
                  paired_mesh_max_axis_spacing_mm=.25, mesh_convergence_checked=False,
                  split_shift_vs_relaxed_percentage_points=100*(candidate['side_share']-ref['side_share']),
                  hits_split_target=bool(.499<=candidate['side_share']<=.501))
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    print(json.dumps(report))


if __name__=='__main__':
    run(Path(sys.argv[1]),Path(sys.argv[2]))
