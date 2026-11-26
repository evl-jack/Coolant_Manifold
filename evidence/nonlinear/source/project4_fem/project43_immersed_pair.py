"""Paired sharp density masks on a shared rectangular Taylor-Hood mesh (Re=50)."""
import argparse
import json
from pathlib import Path

import numpy as np

from .fem import LOCAL_EDGES, TaylorHood
from .junction import strip_triangles
from .project43_explicit import axis
from .project43_leakage import directional_integrals
from .project43_power import port_power


CASES = dict(balanced=(0.814874274531841, 0.), revised=(0.9885282260968125, 0.7528679434757969))


def leakage(fem, fluid, velocity, inlet, throat):
    faces = {}
    for i, tri in enumerate(fem.triangles):
        for slot, (a, b) in enumerate(LOCAL_EDGES):
            v, w = tri[a], tri[b]
            faces.setdefault(tuple(sorted((v, w))), []).append((i,v,w,fem.velocity_dofs[i,3+slot]))
    gross = np.zeros(2)
    length = 0.
    for adj in faces.values():
        if len(adj) != 2 or fluid[adj[0][0]] == fluid[adj[1][0]]:
            continue
        _,v,w,m = next(edge for edge in adj if fluid[edge[0]])
        dx,dy = fem.points[w]-fem.points[v]
        length += np.hypot(dx,dy)
        gross += directional_integrals(*(velocity[[v,m,w]] @ [dy,-dx]))
    np.testing.assert_allclose(length*4,25+2*(4-throat))
    return dict(interface_mm=length*4, fluid_to_solid_fraction=gross[0]/inlet,
                solid_to_fluid_fraction=gross[1]/inlet, net_into_solid_fraction=(gross[0]-gross[1])/inlet)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--spacing',type=float,default=.125)
    parser.add_argument('--alpha',type=float,default=1e6)
    args = parser.parse_args()
    if min(args.alpha,args.spacing)<=0 or not np.all(np.isfinite([args.alpha,args.spacing])):
        parser.error('spacing and alpha must be positive finite numbers')
    args.output.mkdir()
    shoulder=CASES['revised'][1]
    x=axis([0,9.5-shoulder,9.5,10.5,10.5+shoulder,13,14,20],args.spacing)
    y=axis([0,CASES['balanced'][0],CASES['revised'][0],4,6,7],args.spacing)
    xx,yy=np.meshgrid(x,y)
    fem=TaylorHood(np.column_stack((xx.ravel(),yy.ravel()))/4,
                   strip_triangles(np.arange(xx.size).reshape(xx.shape)))
    edges=fem.boundary_edges
    xm,ym=(4*fem.velocity_points[edges[:,2]]).T
    inlet=np.isclose(xm,0)&(ym<4)
    main_port=np.isclose(xm,20)&(ym<4)
    side=np.isclose(ym,7)&(xm>9.5)&(xm<10.5)
    ports=edges[inlet],edges[main_port],edges[side]
    widths=[4*np.linalg.norm(fem.points[p[:,1]]-fem.points[p[:,0]],axis=1).sum() for p in ports]
    np.testing.assert_allclose(widths,[4,4,1])
    def profile(pos):
        xn,yn=pos.T
        return np.column_stack((np.where(np.isclose(xn,0)&(yn<1),6*yn*(1-yn),0),np.zeros_like(yn)))
    fixed,values=fem.velocity_boundary(np.unique(edges[~(main_port|side)]),profile)
    flux=np.array([fem.flux_functional(p,n) for p,n in zip(ports,[(1,0),(1,0),(0,1)])])
    xq,yq=4*fem.quad_points.transpose(2,0,1)
    mu,rho,H=.001002,998.2,.004
    speed=50*mu/(rho*H)
    report=dict(reynolds=50,alpha_kg_m3_s=args.alpha,spacing_mm=args.spacing,
                penetration_mm=1000*np.sqrt(mu/args.alpha),cases={})
    for name,(h,s) in CASES.items():
        gamma=((yq<4)&~((xq>13)&(xq<14)&(yq>h)))|((yq>4)&(abs(xq-10)<.5+s*(yq<6)))
        assert np.all(gamma==gamma[:,:1])
        area=16*np.sum(fem.weights*gamma)
        np.testing.assert_allclose(area,79+h+4*s)
        drag=args.alpha*H**2/mu*(~gamma)
        state=None
        for re in (0.,25.,50.):
            fem.inertia=re
            state=fem.solve(fixed,values,fem.load(),resistance=drag,initial=None if state is None else state.values)
        w=state.values
        q=flux@w
        balance=abs(q[0]-sum(q[1:]))/q[0]
        assert min(q)>0 and balance<1e-9
        vel=w[:2*fem.n_u].reshape(2,-1).T
        power=sum(port_power(fem,ports,w,50.)[0])*mu*speed**2
        result=dict(open_height_mm=h,shoulder_width_mm=s,area_mm2=area,side_share=q[2]/q[0],
                    power_W_per_m=power,mass_imbalance=balance,final_residual=state.residual_history[-1],
                    leakage=leakage(fem,gamma[:,0],vel,q[0],h))
        report['cases'][name]=result
        np.savez_compressed(args.output/(name+'.npz'),state=w,vertices_mm=fem.points*4,
                            triangles=fem.triangles,velocity_dofs=fem.velocity_dofs,fluid=gamma[:,0])
        print(name,result,flush=True)
    a,b=report['cases'].values()
    report.update(split_change_pp=100*(b['side_share']-a['side_share']),power_saving_fraction=1-b['power_W_per_m']/a['power_W_per_m'])
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
