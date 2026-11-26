"""Unoptimized, boundary-fitted 2D laminar T-junction. Run with python -m project4_fem.project43_baseline."""

import argparse
from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np

from .fem import TaylorHood
from .fillet import pressure_weights
from .junction import boundary_data, strip_triangles


@dataclass(frozen=True)
class Geometry:
    length_mm: float = 20.
    height_mm: float = 4.
    branch_center_mm: float = 10.
    branch_width_mm: float = 1.
    branch_length_mm: float = 3.


def fitted_mesh(geometry, spacing_mm):
    g = geometry
    left = g.branch_center_mm - g.branch_width_mm / 2
    right = g.branch_center_mm + g.branch_width_mm / 2
    if not all(np.isfinite(v) and v > 0 for v in (*asdict(g).values(), spacing_mm)):
        raise ValueError("Geometry and mesh spacing must be finite and positive")
    if not 0 < left < right < g.length_mm:
        raise ValueError("Branch mouth must lie strictly between main port planes")

    def segment(start, end):
        intervals = max(1, int(np.ceil((end - start) / spacing_mm)))
        return np.linspace(start, end, intervals + 1)

    upstream = segment(0., left)
    mouth = segment(left, right)
    main_x = np.r_[upstream[:-1], mouth, segment(right, g.length_mm)[1:]]
    main_y = segment(0., g.height_mm)
    branch_y = segment(g.height_mm, g.height_mm + g.branch_length_mm)[1:]

    def grid(x, y):
        xx, yy = np.meshgrid(x, y)
        return np.column_stack((xx.ravel(), yy.ravel()))

    main = grid(main_x, main_y)
    branch = grid(mouth, branch_y)
    main_ids = np.arange(len(main)).reshape(len(main_y), len(main_x))
    start = len(upstream) - 1
    branch_ids = np.vstack((main_ids[-1, start:start+len(mouth)],
                            np.arange(len(main), len(main)+len(branch)).reshape(len(branch_y), len(mouth))))
    points = np.vstack((main, branch)) / g.height_mm
    return points, np.vstack((strip_triangles(main_ids), strip_triangles(branch_ids)))


def baseline(geometry, spacing_mm, reynolds):
    """Re is based on inlet height and bulk velocity. All fluxes are per unit depth."""
    if not np.isfinite(reynolds) or reynolds <= 0:
        raise ValueError("Reynolds number must be finite and positive")
    fem = TaylorHood(*fitted_mesh(geometry, spacing_mm))
    fixed, prescribed, ports, flux = boundary_data(fem)
    port_widths = np.array([
        np.linalg.norm(fem.points[edges[:, 1]] - fem.points[edges[:, 0]], axis=1).sum()
        for edges in ports]) * geometry.height_mm
    area_mm2 = fem.weights.sum() * geometry.height_mm**2
    expected_area = geometry.length_mm*geometry.height_mm + geometry.branch_width_mm*geometry.branch_length_mm
    np.testing.assert_allclose(area_mm2, expected_area, rtol=1e-12)
    np.testing.assert_allclose(port_widths, [geometry.height_mm, geometry.height_mm, geometry.branch_width_mm],
                               rtol=1e-12, atol=0)

    load = fem.load()
    state = None
    continuation = []
    for inertia in (0., reynolds/2, reynolds):
        fem.inertia = inertia
        state = fem.solve(fixed, prescribed, load, initial=None if state is None else state.values)
        continuation.append(dict(inertia=inertia, residuals=state.residual_history))
    q_in, q_main, q_side = flux @ state.values
    imbalance = abs(q_main + q_side - q_in) / q_in
    if not (abs(q_in - 1) < 1e-10 and imbalance < 1e-10 and min(q_main, q_side) > 0):
        raise RuntimeError(f"Baseline failed signed flux checks: {q_in}, {q_main}, {q_side}")
    p_in, p_out = pressure_weights(fem, ports)
    loss = p_in @ state.values - (state.values @ (p_out @ state.values)) / (q_main + q_side)
    rho, mu = 998.2, 0.001002
    height_m = geometry.height_mm / 1000
    speed = reynolds*mu / (rho*height_m)
    flow_scale = speed*height_m
    pressure_scale = mu*speed / height_m
    report = dict(
        case="Unoptimized, clear-fluid, fitted-wall laminar Navier-Stokes",
        geometry_mm=asdict(geometry), fluid_area_mm2=float(area_mm2),
        port_widths_mm=dict(zip(("inlet", "main", "side"), port_widths.tolist())),
        mesh_max_axis_spacing_mm=spacing_mm, triangles=len(fem.triangles), dofs=fem.size,
        fluid_density_kg_m3=rho, dynamic_viscosity_Pa_s=mu, reynolds_height=reynolds,
        bulk_inlet_speed_m_s=speed, reference_height_m=height_m, pressure_scale_Pa=pressure_scale,
        boundary_conditions="Fully developed parabolic inlet, no-slip walls, zero symmetric traction at both exits",
        flux_per_depth_m2_s=dict(inlet=q_in*flow_scale, main=q_main*flow_scale, side=q_side*flow_scale),
        side_share=q_side/q_in, side_to_main=q_side/q_main, relative_mass_imbalance=imbalance,
        pressure_loss_Pa=loss*pressure_scale,
        pressure_loss_definition="Inlet area-mean static pressure minus signed flux-weighted static pressure over both outlets",
        continuation=continuation)
    return report, fem, state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, value in asdict(Geometry()).items():
        parser.add_argument("--"+name.replace("_", "-"), type=float, default=value)
    parser.add_argument("--spacing-mm", type=float, default=0.25)
    parser.add_argument("--reynolds", type=float, default=50.)
    parser.add_argument("--output", type=Path, help="New directory for report.json and field.npz; defaults to stdout only")
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error("Choose a new output directory to preserve previous results")
    g = Geometry(**{key: getattr(args, key) for key in asdict(Geometry())})
    report, fem, state = baseline(g, args.spacing_mm, args.reynolds)
    formatted = json.dumps(report, indent=2, allow_nan=False)
    if args.output:
        args.output.mkdir(parents=True)
        np.savez_compressed(args.output/"field.npz", vertices_mm=fem.points*g.height_mm,
                            triangles=fem.triangles, velocity_nodes_mm=fem.velocity_points*g.height_mm,
                            velocity_dofs=fem.velocity_dofs, velocity_m_s=state.values[:2*fem.n_u].reshape(2, -1).T
                            * report["bulk_inlet_speed_m_s"], pressure_Pa=state.values[2*fem.n_u:]
                            * report["pressure_scale_Pa"])
        (args.output/"report.json").write_text(formatted+"\n", encoding="utf-8")
    print(formatted)


if __name__ == "__main__":
    main()
