"""Frozen binary Brinkman T-junction on a rectangular background, using viscous pressure units.

Run from the package parent with python -m project4_fem.project43_brinkman --help.
The original bottom wall is fitted. The main roof and the two branch sides are
immersed interfaces. Solid parts of the exterior boundary are closed.
"""

import argparse
from dataclasses import asdict
import json
from pathlib import Path

import numpy as np

from .fem import TaylorHood
from .fillet import pressure_weights
from .junction import strip_triangles
from .project43_baseline import Geometry, fitted_mesh


def create_problem(g, spacing_mm):
    fitted, _ = fitted_mesh(g, spacing_mm)
    x, y = np.unique(fitted[:, 0]), np.unique(fitted[:, 1])
    xx, yy = np.meshgrid(x, y)
    nodes = np.arange(xx.size).reshape(xx.shape)
    fem = TaylorHood(np.column_stack((xx.ravel(), yy.ravel())), strip_triangles(nodes))
    left = (g.branch_center_mm - g.branch_width_mm / 2) / g.height_mm
    right = (g.branch_center_mm + g.branch_width_mm / 2) / g.height_mm
    xq, yq = fem.quad_points.transpose(2, 0, 1)
    gamma = ((yq < 1.) | ((xq > left) & (xq < right))).astype(float)
    if np.any(gamma != gamma[:, :1]):
        raise ValueError("Frozen interface must follow triangle edges")

    edges = fem.boundary_edges
    xm, ym = fem.velocity_points[edges[:, 2]].T
    inlet = np.isclose(xm, x[0]) & (ym < 1.)
    main = np.isclose(xm, x[-1]) & (ym < 1.)
    side = np.isclose(ym, y[-1]) & (left < xm) & (xm < right)
    ports = edges[inlet], edges[main], edges[side]
    normals = ((1, 0), (1, 0), (0, 1))
    flux = np.array([fem.flux_functional(port, normal) for port, normal in zip(ports, normals)])

    def profile(points):
        xnode, eta = points.T
        u = np.where(np.isclose(xnode, 0.) & (eta < 1.), 6 * eta * (1 - eta), 0.)
        return np.column_stack((u, np.zeros_like(u)))

    fixed, values = fem.velocity_boundary(np.unique(edges[~(main | side)]), profile)
    widths = [np.linalg.norm(fem.points[p[:, 1]] - fem.points[p[:, 0]], axis=1).sum()
              * g.height_mm for p in ports]
    np.testing.assert_allclose(widths, [g.height_mm, g.height_mm, g.branch_width_mm], rtol=1e-12)
    return fem, gamma, fixed, values, ports, flux, widths


def run(geometry, spacing_mm, reynolds, alpha_max, check_gradient=True):
    if not all(np.isfinite(v) and v > 0 for v in (reynolds, alpha_max)):
        raise ValueError("Reynolds number and dimensional penalty must be finite and positive")
    fem, gamma, fixed, values, ports, flux, widths = create_problem(geometry, spacing_mm)
    rho, mu = 998.2, 0.001002
    height = geometry.height_mm * 1e-3
    bulk = reynolds * mu / (rho * height)
    pressure_scale = mu * bulk / height
    alpha_star = alpha_max * height**2 / mu
    resistance = alpha_star * (1 - gamma)**2
    areas = [np.sum(fem.weights * mask) * geometry.height_mm**2 for mask in (gamma, 1 - gamma)]
    expected_fluid = geometry.length_mm * geometry.height_mm + geometry.branch_width_mm * geometry.branch_length_mm
    expected_solid = (geometry.length_mm - geometry.branch_width_mm) * geometry.branch_length_mm
    np.testing.assert_allclose(areas, [expected_fluid, expected_solid], rtol=1e-12)

    load = fem.load()
    p_in, p_out = pressure_weights(fem, ports)
    state = None
    continuation = []
    for inertia in (0., reynolds / 2, reynolds):
        fem.inertia = inertia
        state = fem.solve(fixed, values, load, resistance=resistance,
                          initial=state.values if state is not None else None)
        continuation.append(dict(inertia=inertia, residuals=state.residual_history))

    def responses(w):
        inlet, main, side = flux @ w
        imbalance = abs(main + side - inlet) / inlet
        if not (abs(inlet - 1) < 1e-10 and imbalance < 1e-10 and min(main, side) > 0):
            raise RuntimeError(f"Invalid inlet flux or outlet balance: {inlet}, {main}, {side}")
        return np.array([side / inlet, pressure_scale * (p_in @ w - w @ (p_out @ w) / (main + side))])

    w = state.values
    side_share, loss_pa = responses(w)
    inlet, main, side = flux @ w
    outflow = main + side
    partials = np.array([(flux[2] - side_share * flux[0]) / inlet,
                         pressure_scale * (p_in - (p_out + p_out.T) @ w / outflow
                                           + (w @ (p_out @ w)) * (flux[1] + flux[2]) / outflow**2)])
    # For s=ln(alpha_max / reference), d(alpha_star * (1-gamma)^2)/ds = resistance.
    adjoint = fem.drag_adjoint(state, partials, resistance[:, :, None])[:, 0]
    sensitivity = dict(parameter="ln(alpha_max / reference) at fixed fluid fraction and mesh",
                       quantities=["side_share_fraction", "pressure_loss_Pa"], adjoint=adjoint.tolist())
    if check_gradient:
        step = 1e-3
        perturbed = [responses(fem.solve(fixed, values, load, resistance=np.exp(k * step) * resistance,
                                         initial=w).values) for k in (-1, 1)]
        centered = (perturbed[1] - perturbed[0]) / (2 * step)
        sensitivity.update(log_step=step, centered=centered.tolist(),
                           absolute_difference=np.abs(adjoint - centered).tolist())
        np.testing.assert_allclose(adjoint, centered, rtol=2e-4, atol=1e-10)

    u, _ = fem.velocity_fields(w)
    solid_weights = (1 - gamma) * fem.weights
    rms = np.sqrt(np.sum(solid_weights * np.sum(u*u, axis=-1)) / solid_weights.sum())
    penetration_mm = 1000 * np.sqrt(mu / alpha_max)
    report = dict(case="Frozen binary T-junction, steady laminar Navier-Stokes-Brinkman",
                  geometry_mm=asdict(geometry), representation=__doc__.splitlines()[-2:],
                  gamma_definition="1 fluid, 0 solid, sampled at quadrature; alpha=alpha_max*(1-gamma)^2",
                  fluid_area_mm2=areas[0], solid_area_mm2=areas[1],
                  port_widths_mm=dict(zip(("inlet", "main", "side"), widths)),
                  mesh_max_axis_spacing_mm=spacing_mm, triangles=len(fem.triangles), dofs=fem.size,
                  fluid_density_kg_m3=rho, dynamic_viscosity_Pa_s=mu, reynolds_height=reynolds,
                  bulk_inlet_speed_m_s=bulk, reference_height_m=height, pressure_scale_Pa=pressure_scale,
                  boundary_conditions="Parabolic inlet, closed exterior wall segments, zero symmetric traction on the two outlet apertures",
                  alpha_max_kg_m3_s=alpha_max, alpha_star=alpha_star, brinkman_length_mm=penetration_mm,
                  brinkman_length_over_mesh_spacing=penetration_mm / spacing_mm,
                  flux_per_depth_m2_s=dict(zip(("inlet", "main", "side"), (flux @ w) * bulk * height)),
                  side_share=side_share, side_to_main=side/main, relative_mass_imbalance=abs(outflow-inlet)/inlet,
                  solid_rms_speed_over_inlet_bulk=rms, solid_rms_speed_m_s=rms * bulk,
                  pressure_loss_Pa=loss_pa,
                  pressure_loss_definition="Inlet area-mean static pressure minus signed flux-weighted static pressure over both outlets",
                  sensitivity=sensitivity, continuation=continuation)
    return report, fem, state, gamma


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key, default in asdict(Geometry()).items():
        parser.add_argument("--" + key.replace("_", "-"), type=float, default=default)
    parser.add_argument("--spacing-mm", type=float, default=0.25)
    parser.add_argument("--reynolds", type=float, default=50.)
    parser.add_argument("--alpha-max", type=float, default=1e6, help="Dynamic penalty in kg/(m³ s)")
    parser.add_argument("--skip-gradient-check", action="store_true", help="Skip the two finite-difference re-solves")
    parser.add_argument("--output", type=Path, help="Fresh directory for report.json and field.npz")
    args = parser.parse_args()
    if args.output and args.output.exists():
        parser.error("Choose a fresh output directory")
    g = Geometry(**{key: getattr(args, key) for key in asdict(Geometry())})
    report, fem, state, gamma = run(g, args.spacing_mm, args.reynolds, args.alpha_max, not args.skip_gradient_check)
    formatted = json.dumps(report, indent=2, allow_nan=False)
    if args.output:
        args.output.mkdir(parents=True)
        np.savez_compressed(args.output / "field.npz", vertices_mm=fem.points*g.height_mm,
                            triangles=fem.triangles, velocity_nodes_mm=fem.velocity_points*g.height_mm,
                            velocity_dofs=fem.velocity_dofs, fluid_fraction_per_triangle=gamma[:, 0],
                            velocity_m_s=state.values[:2*fem.n_u].reshape(2, -1).T*report["bulk_inlet_speed_m_s"],
                            pressure_Pa=state.values[2*fem.n_u:]*report["pressure_scale_Pa"])
        (args.output / "report.json").write_text(formatted + "\n", encoding="utf-8")
    print(formatted)


if __name__ == "__main__":
    main()
