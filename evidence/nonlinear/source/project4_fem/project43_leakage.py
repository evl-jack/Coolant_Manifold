"""Integrate positive and negative P2 normal flow on fitted density interfaces."""
import json
from pathlib import Path
import sys
import numpy as np
from .fem import LOCAL_EDGES


def directional_integrals(v0, vm, v1):
    a = 2 * (v1 + v0 - 2*vm)
    b = v1 - v0 - a
    roots = np.roots([a, b, v0])
    cuts = sorted([0., 1.] + [float(r.real) for r in roots if abs(r.imag) < 1.e-13 and 0 < r.real < 1])
    primitive = lambda t: t*(v0 + t*(b/2 + t*a/3))
    increments = np.diff([primitive(t) for t in cuts])
    return np.maximum(increments, 0).sum(), -np.minimum(increments, 0).sum()


def run(directory):
    report = json.loads((directory/'report.json').read_text())
    with np.load(directory/'field.npz') as field:
        triangles, points = field['triangles'], field['vertices_mm']*1.e-3
        dofs, velocity = field['velocity_dofs'], field['velocity_m_s']
        fluid = field['fluid_fraction_per_triangle']
    assert np.all((fluid == 0) | (fluid == 1))
    edges = {}
    for i, tri in enumerate(triangles):
        for slot, (a, b) in enumerate(LOCAL_EDGES):
            v, w = tri[a], tri[b]
            edges.setdefault(tuple(sorted((v, w))), []).append((i, v, w, dofs[i, 3+slot]))
    totals = np.zeros(2)
    length = 0.
    for adjacent in edges.values():
        if len(adjacent) != 2 or fluid[adjacent[0][0]] == fluid[adjacent[1][0]]:
            continue
        i, v, w, m = next(edge for edge in adjacent if fluid[edge[0]] == 1)
        dx, dy = points[w] - points[v]
        length += np.hypot(dx, dy)
        # Counter-clockwise cell orientation gives the outward length-weighted normal.
        normal_flux = velocity[[v, m, w]] @ [dy, -dx]
        totals += directional_integrals(*normal_flux)
    g = report['geometry_mm']
    np.testing.assert_allclose(length*1000, g['length_mm']-g['branch_width_mm']+2*g['branch_length_mm'])
    inlet = report['flux_per_depth_m2_s']['inlet']
    assert inlet > 0
    return dict(source=str(directory), interface_length_mm=length*1000,
                fluid_to_solid_fraction=totals[0]/inlet, solid_to_fluid_fraction=totals[1]/inlet,
                net_into_solid_fraction=(totals[0]-totals[1])/inlet, inlet_flux_m2_s=inlet)


if __name__ == "__main__":
    print(json.dumps(run(Path(sys.argv[1])), indent=2, allow_nan=False))
