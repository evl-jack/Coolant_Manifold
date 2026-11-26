"""Incompressible SST-2003m flat-plate screen against TMR VERIF/2DZP.

Reference: https://tmbwg.github.io/turbmodels/flatplate_sst.html
This custom coarse triangular grid uses a uniform far-field Dirichlet trace
and a natural traction outlet. TMR's quoted comparator is compressible SST-Vm
on its 273x193 grid. A single discrepancy therefore measures aggregate error.
"""
import argparse
import json
from pathlib import Path

import numpy as np

from .fem import TaylorHood, quadratic_basis
from .next_rans_adjoint.sst import FittedRANS
from .project43_sst import FrozenSST, pseudo_transient


class FlatPlate(FittedRANS):
    reaction_coefficients = FrozenSST.reaction_coefficients

    def element(self, coefficients, geometry):
        residual = super().element(coefficients, geometry)
        rates = self.reaction_coefficients(coefficients, geometry)
        for block, rate in zip((slice(15, 18), slice(18, 21)), rates[:2]):
            z = coefficients[:, block]
            difference = z[:, None, :] - (z @ self.p.T)[:, :, None]
            residual[:, block] += np.einsum("eq,qa,eqa->ea", geometry[3]*rate, self.p, difference)
        return residual


def run(viscous_wall_start=False):
    reynolds = 5.e6
    k_inf, omega_inf = 1.125/reynolds, 125.
    xs = np.r_[np.linspace(-1/3, 0., 7)[:-1], 2.*np.linspace(0., 1., 41)**1.5]
    ys = np.r_[0., np.geomspace(4.e-6, 1., 43)]
    vertices = np.stack(np.meshgrid(xs, ys), axis=-1).reshape(-1, 2)
    ids = np.arange(len(vertices)).reshape(len(ys), len(xs))
    a, b, c, d = (v.ravel() for v in (ids[:-1,:-1], ids[:-1,1:], ids[1:,1:], ids[1:,:-1]))
    tris = np.vstack((np.c_[a,b,c], np.c_[a,c,d]))
    mesh = TaylorHood(vertices, tris)
    edge = mesh.boundary_edges[:, :2]
    endpoints = vertices[edge]
    plate = edge[np.all(endpoints[:,:,1] == 0., axis=1) & np.all(endpoints[:,:,0] >= 0., axis=1)]
    m = FlatPlate(vertices, tris, reynolds, plate)
    m.solid = np.zeros(len(tris), dtype=bool)
    xv, yv = mesh.velocity_points.T
    x, y = vertices.T
    wall = (y == 0.) & (x >= 0.)
    outer = (y == 1.) | (x == xs[0])
    wall_v = (yv == 0.) & (xv >= 0.)
    outer_v = (yv == 1.) | (xv == xs[0])
    idx = np.flatnonzero
    fixed = np.r_[idx(wall_v | outer_v), mesh.n_u+idx((yv == 0.) | outer_v),
                  mesh.size+idx(wall | outer), mesh.size+m.nb+idx(wall | outer)]
    delta = .37*np.maximum(x, .0001)**.8/reynolds**.2
    delta_v = .37*np.maximum(xv, .0001)**.8/reynolds**.2
    ux = np.where(xv < 0., 1., np.minimum(1., (yv/delta_v)**(1./7.)))
    if viscous_wall_start:
        # Bound the initial wall slope using a local flat-plate friction estimate.
        local_re = reynolds * np.maximum(xv, .0001)
        initial_cf = .0592 * local_re**(-.2)
        ux = np.where(xv < 0., 1., np.minimum(ux, .5 * reynolds * initial_cf * yv))
    ux[outer_v], ux[wall_v] = 1., 0.
    k = k_inf + .003*(1.-np.exp(-y/2.e-5))**2*np.exp(-y/delta)*(x >= 0.)
    omega = omega_inf + 6./(.075*reynolds*np.maximum(y, ys[1])**2)*(x >= 0.)
    k[outer], omega[outer] = k_inf, omega_inf
    k[wall], omega[wall] = 0., 60./(.075*reynolds*ys[1]**2)
    z = np.r_[ux, np.zeros(mesh.n_u+m.nb), k, omega]
    free = np.setdiff1d(np.arange(m.size), fixed)
    r0 = m.residual(z)
    scale = np.ones(m.size)
    for block in m.blocks:
        indices = free[(free >= block.start) & (free < block.stop)]
        scale[block] = max(1., np.max(np.abs(r0[indices])))
    z, history, failure = pseudo_transient(m, z, free, scale, max_steps=100)
    # Evaluate physical wall gradient using the P2 element adjacent to the plate.
    ix = np.searchsorted(xs, .97008)-1
    xi = (.97008-xs[ix])/(xs[ix+1]-xs[ix])
    cell = ix
    _, grad_ref = quadratic_basis(np.array([[1.-xi, xi, 0.]]))
    jac = (vertices[tris[cell,1:]]-vertices[tris[cell,0]]).T
    grad = grad_ref[0] @ np.linalg.inv(jac)
    du = z[mesh.velocity_dofs[cell]] @ grad
    dv = z[mesh.n_u+mesh.velocity_dofs[cell]] @ grad
    cf = 2.*(du[1]+dv[0])/reynolds
    r = m.residual(z)
    return dict(converged=failure is None, failure=failure, history=history,
                reynolds_per_length=reynolds, triangle_count=len(tris), first_dy=ys[1],
                variant="SST-2003m, incompressible, lumped destruction",
                viscous_wall_start=viscous_wall_start,
                max_scaled_residual=float(np.max(np.abs(r[free]/scale[free]))),
                probe_x=.97008, computed_cf=float(cf) if failure is None else None,
                reference_cf=.00268299, reference="TMR SST-Vm CFL3D 273x193, Mach 0.2",
                relative_cf_difference=float(cf/.00268299-1.) if failure is None else None)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--viscous-wall-start", action="store_true")
    args = parser.parse_args()
    result = run(viscous_wall_start=args.viscous_wall_start)
    with args.output.open("x", encoding="utf8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
    raise SystemExit(0 if result["converged"] else 2)
