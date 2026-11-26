"""Boundary-fitted diagnostic for the reference and a constrained feasibility probe."""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

from .fem import TaylorHood
from .fillet import BASELINE, FilletMap, pressure_weights
from .junction import boundary_data, junction_mesh


def evaluate(geometry, design, refine):
    mapping = FilletMap(geometry)
    geometry_offset = lambda y: mapping.offset(y, design)
    fem = TaylorHood(*junction_mesh(geometry, refine, geometry_offset))
    fixed, values, ports, flux = boundary_data(fem)
    p_in, p_out = pressure_weights(fem, ports)
    load = fem.load()
    viscous = fem.viscous_pressure[:2*fem.n_u, :2*fem.n_u]
    evaluations = []
    initial = None
    for re in (0., 25., 50.):
        fem.inertia = re
        state = fem.solve(fixed, values, load, initial=initial)
        initial = state.values
        q, qm, qs = flux@initial
        u = initial[:2*fem.n_u]
        evaluations.append(dict(reynolds=re, side_share=float(qs/q), side_to_main=float(qs/qm),
                                pressure_loss=float(p_in@initial-initial@(p_out@initial)/(qm+qs)),
                                dissipation=float(u@(viscous@u)), mass_imbalance=float(abs(q-qm-qs)/q),
                                residual=state.residual_history[-1]))

    y = np.linspace(mapping.root, mapping.root+mapping.radius, 1001)
    width = 2*(mapping.stem+mapping.parameters(y)[1]+geometry_offset(y))
    area = 2*quad(lambda y: mapping.offset(np.asarray(y)) - geometry_offset(np.asarray(y)),
                  mapping.root, mapping.root+mapping.radius, epsabs=1e-8, limit=100)[0]
    return dict(triangles=len(fem.triangles), minimum_width_mm=float(width.min()),
                geometric_solid_area_mm2=area, responses=evaluations)


def run(study, output):
    limits = json.loads((study/"limits.json").read_text())
    report = json.loads((study/"search-re-50.json").read_text())
    geometry = (1., .625, .5, 1.5)
    entries = {}
    for name, design in (("baseline", BASELINE), ("pressure-and-area-feasibility-probe", report["probe"]["design"])):
        meshes = [evaluate(geometry, np.array(design), level) for level in (1, 2)]
        entries[name] = dict(design=list(design), meshes=meshes,
                             side_share_grid_change_pp=100*(meshes[1]["responses"][-1]["side_share"]
                                                           -meshes[0]["responses"][-1]["side_share"]))
        print(f"{name}: explicit fine-grid Re 50 side share {meshes[-1]['responses'][-1]['side_share']:.4%}", flush=True)
    data = dict(target=limits["target"], tolerance=limits["side_tolerance"], geometry_mm=list(geometry), **entries)
    output.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.study, args.output)
