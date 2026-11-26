"""Boundary-fitted T-junction reference for the reduced Project 4 model."""

import argparse
import json
from pathlib import Path
from time import perf_counter

import numpy as np

from .fem import TaylorHood


def strip_triangles(nodes):
    a, b, c, d = (n.ravel() for n in (nodes[:-1, :-1], nodes[:-1, 1:], nodes[1:, 1:], nodes[1:, :-1]))
    return np.vstack((np.column_stack((a, b, c)), np.column_stack((a, c, d))))


def junction_mesh(controls, refine, fillet_offset=None, transition_rows=6):
    """Join a rectangular main passage to a fitted strip in the branch."""
    main, outlet, stem, fillet = controls
    bottom, root, center, guide, exit_y = 4-main, 4+main, 10., 8., 11.
    length_scale = 2*main
    columns = (8 if fillet_offset is None else 16)*refine
    s = np.linspace(-1., 1., columns+1)
    if fillet_offset is not None:
        s = np.sin(.5*np.pi*s)
    mouth = center + s*(stem+fillet)
    x = np.r_[np.linspace(-1., mouth[0], 13*refine+1)[:-1], mouth,
              np.linspace(mouth[-1], 21., 13*refine+1)[1:]]
    y = np.linspace(bottom, root, 5*refine+1)
    main_points = np.column_stack([a.ravel() for a in np.meshgrid(x, y)])
    main_nodes = np.arange(len(main_points)).reshape(len(y), len(x))

    relative = np.linspace(0., 1., (8 if fillet_offset is None else 38)*refine+1)
    if fillet_offset is not None:
        relative = .4*relative + .6*(1-np.cos(.5*np.pi*relative))
    ys = np.r_[root+fillet*relative[1:],
               np.linspace(root+fillet, guide, transition_rows*refine+1)[1:],
               np.linspace(guide, exit_y, 9*refine+1)[1:]]
    h = np.minimum(ys-root, fillet)
    t = np.clip((ys-root-fillet)/(guide-root-fillet), 0, 1)
    width = stem + fillet - np.sqrt(h*(2*fillet-h)) + (outlet-stem)*t*t*(3-2*t)
    if fillet_offset is not None:
        width += fillet_offset(ys)
    branch_points = np.stack((center + width[:, None]*s, np.broadcast_to(ys[:, None], (len(ys), len(s)))), axis=-1).reshape(-1, 2)
    branch_nodes = np.vstack((main_nodes[-1, 13*refine:13*refine+len(s)],
                              np.arange(len(main_points), len(main_points)+len(branch_points)).reshape(len(ys), len(s))))
    points = (np.vstack((main_points, branch_points))-[-1., bottom])/length_scale
    triangles = np.vstack((strip_triangles(main_nodes), strip_triangles(branch_nodes)))
    return points, triangles


def boundary_data(fem):
    edges, points = fem.boundary_edges, fem.velocity_points
    mid = points[edges[:, 2]]
    inlet = np.isclose(mid[:, 0], 0.)
    main = np.isclose(mid[:, 0], points[:, 0].max())
    side = np.isclose(mid[:, 1], points[:, 1].max())
    no_slip_or_inlet = np.unique(edges[~(main | side)])

    def boundary_velocity(position):
        x, y = position.T
        u = np.where(np.isclose(x, 0.), 6*y*(1-y), 0.)
        return np.column_stack((u, np.zeros_like(u)))

    fixed, values = fem.velocity_boundary(no_slip_or_inlet, boundary_velocity)
    ports = edges[inlet], edges[main], edges[side]
    flux_rows = [fem.flux_functional(port, normal) for port, normal in zip(ports, ((1, 0), (1, 0), (0, 1)))]
    return fixed, values, ports, np.array(flux_rows)


def reference_case(controls, refinement, reynolds, output):
    start = perf_counter()
    fem = TaylorHood(*junction_mesh(controls, refinement))
    fixed, values, ports, flux_rows = boundary_data(fem)
    points = fem.velocity_points
    load = fem.load()
    results, initial = [], None
    for inertia in (0., reynolds/2, reynolds):
        fem.inertia = inertia
        state = fem.solve(fixed, values, load, initial=initial)
        initial = state.values
        q_in, q_main, q_side = (float(row @ state.values) for row in flux_rows)
        pressure = state.values[2*fem.n_u:]
        inlet_edges = ports[0][:, :2]
        inlet_widths = np.linalg.norm(points[inlet_edges[:, 1]]-points[inlet_edges[:, 0]], axis=1)
        inlet_pressure = np.average(pressure[inlet_edges].mean(axis=1), weights=inlet_widths)
        results.append(dict(reynolds=inertia, inlet=q_in, main=q_main, side=q_side,
                            side_share=q_side/q_in, side_to_main=q_side/q_main,
                            mass_imbalance=abs(q_main+q_side-q_in)/q_in,
                            mean_inlet_pressure=float(inlet_pressure), newton_residuals=state.residual_history))
    np.savez_compressed(output / f"reference-grid-{refinement}.npz", points=fem.points, triangles=fem.triangles,
                        velocity_points=points, velocity_dofs=fem.velocity_dofs, state=state.values)
    return dict(refinement=refinement, triangles=len(fem.triangles), state_dofs=fem.size,
                wall_seconds=perf_counter()-start, solutions=results)


def run(checkpoint, reynolds, output):
    source = json.loads(checkpoint.read_text(encoding="utf-8"))["final"]
    names = ("main_radius_mm", "outlet_radius_mm", "stem_radius_mm", "fillet_radius_mm")
    controls = [source[name] for name in names]
    output.mkdir(parents=True, exist_ok=True)
    grids = [reference_case(controls, i, reynolds, output) for i in (1, 2)]
    report = dict(geometry_mm=dict(zip(names, controls)), inlet="parabolic profile, unit flow per depth",
                  outlets="natural zero traction", reference="height and mean inlet velocity",
                  models=grids, side_share_grid_change_pp=100*(grids[1]["solutions"][-1]["side_share"]
                                                            -grids[0]["solutions"][-1]["side_share"]))
    (output / "junction-reference.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--reynolds", type=float, default=50.)
    parser.add_argument("--output", type=Path, default=Path("project4_fem/output"))
    args = parser.parse_args()
    if not np.isfinite(args.reynolds) or args.reynolds <= 0:
        parser.error("Reynolds number must be finite and positive")
    run(args.checkpoint, args.reynolds, args.output)
