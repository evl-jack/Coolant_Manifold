"""Archive an independent hard-wall reconstruction of a topology/baffle endpoint."""

import json
import os
from pathlib import Path

import numpy as np

from ..fem import quadratic_basis
from ..project43_leakage import directional_integrals
from .domain import REFERENCE_HEIGHT_MM
from .evaluator import HEIGHT_M, VELOCITY_SCALE, VISCOSITY, response_values
from .full_density_continuation import (
    FITTED_POWER_RELATIVE_CHANGE, FITTED_SPLIT_CHANGE_PP, FLOW_BAND, MAX_AREA_MM2,
)
from .reconstruction import _PORTS_MM, Reconstruction, reconstruct, solve_fitted, write_reconstruction


def _exact_ports(geometry, contour_segments_mm):
    from shapely.geometry import LineString

    if geometry.is_empty or not geometry.is_valid or geometry.geom_type not in ("Polygon", "MultiPolygon"):
        raise ValueError("The exact strip and baffle operations did not yield valid polygonal fluid")
    polygons = (geometry,) if geometry.geom_type == "Polygon" else tuple(geometry.geoms)
    polygons = tuple(sorted(polygons, key=lambda polygon: (*polygon.bounds, polygon.area)))
    port_components = []
    for port in _PORTS_MM:
        line = LineString(port)
        match = [i for i, polygon in enumerate(polygons) if polygon.boundary.covers(line)]
        if len(match) != 1:
            raise ValueError(f"The exact strip or baffle altered port {port}")
        port_components.append(match[0])
    return Reconstruction(geometry, polygons, tuple(port_components), contour_segments_mm)


def _reconstruction_pair(grid, projected, baffle_height_mm):
    from shapely.geometry import box

    if not np.isfinite(baffle_height_mm) or not 0.0 <= baffle_height_mm <= 3.75:
        raise ValueError("Baffle height must lie in [0, 3.75] mm")
    raw_shape = reconstruct(grid, projected)
    strip = box(13.0, 0.0, 14.0, 4.0)
    pre = _exact_ports(raw_shape.geometry.union(strip), raw_shape.contour_segments_mm)
    if baffle_height_mm == 0.0:
        return pre, pre
    post = pre.geometry.difference(box(13.0, 0.0, 14.0, baffle_height_mm))
    after = _exact_ports(post, raw_shape.contour_segments_mm)
    np.testing.assert_allclose(pre.area_mm2 - after.area_mm2, baffle_height_mm, rtol=1e-12, atol=1e-12)
    return pre, after


def reconstruct_hybrid(grid, projected, baffle_height_mm):
    """Return the gamma=0.5 polygon with the exact strip and hard baffle."""
    return _reconstruction_pair(grid, projected, baffle_height_mm)[1]


def _fitted_checks(geometry, spacings):
    if not spacings:
        return {"attempted": False, "reason": "No fitted mesh sequence was requested"}, None
    if geometry.port_component is None:
        return {"attempted": False, "reason": "The hard-wall polygon disconnects the ports"}, None
    if geometry.area_mm2 > MAX_AREA_MM2:
        return {"attempted": False, "reason": "The hard-wall polygon exceeds the area cap"}, None

    levels = []
    finest = None
    for spacing in spacings:
        try:
            result, fem, state = solve_fitted(geometry, spacing)
        except Exception as error:
            return {"attempted": True, "completed": False, "failed_spacing_mm": spacing,
                    "error": f"{type(error).__name__}: {error}", "levels": levels}, finest
        levels.append(result)
        finest = fem, state
    changes = []
    for coarse, fine in zip(levels[:-1], levels[1:]):
        split = abs(fine["side_share"] - coarse["side_share"])
        power = fine["mechanical_power_W_per_m"]
        relative_power = abs(fine["mechanical_power_W_per_m"] - coarse["mechanical_power_W_per_m"]) / power if power > 0 else None
        changes.append({"split_change_pp": 100 * split, "relative_power_change": relative_power,
                        "passed": 100 * split <= FITTED_SPLIT_CHANGE_PP
                        and relative_power is not None and relative_power <= FITTED_POWER_RELATIVE_CHANGE})
    uncertainty = max((change["split_change_pp"] / 100 for change in changes), default=None)
    side_share = levels[-1]["side_share"]
    interval = None if uncertainty is None else (side_share - uncertainty, side_share + uncertainty)
    verified = (len(levels) >= 3 and all(level["state_passed"] for level in levels)
                and all(change["passed"] for change in changes)
                and FLOW_BAND[0] <= interval[0] <= interval[1] <= FLOW_BAND[1])
    return {"attempted": True, "completed": True, "levels": levels, "successive_changes": changes,
            "physical_flow_band": FLOW_BAND, "conservative_split_uncertainty": uncertainty,
            "conservative_fitted_split_interval": interval, "fitted_diagnostic_passed": bool(verified)}, finest


def _state_fields(fem, state):
    values = np.asarray(state.values)
    if values.shape != (fem.size,) or not np.isfinite(values).all():
        raise ValueError("Expected a finite state vector on the supplied FEM mesh")
    velocity = VELOCITY_SCALE * values[:2 * fem.n_u].reshape((2, fem.n_u)).T
    pressure = VISCOSITY * VELOCITY_SCALE / HEIGHT_M * values[2 * fem.n_u:]
    return velocity, pressure


def _write_state(output, stem, fem, state):
    velocity, pressure = _state_fields(fem, state)
    with (output / f"{stem}.npz").open("xb") as stream:
        np.savez_compressed(stream, points_mm=fem.points * REFERENCE_HEIGHT_MM,
                            triangles=fem.triangles, quadratic_points_mm=fem.velocity_points * REFERENCE_HEIGHT_MM,
                            velocity_dofs=fem.velocity_dofs, velocity_m_s=velocity, pressure_pa=pressure)
    speed = np.linalg.norm(velocity, axis=1)
    return {"archive": f"{stem}.npz", "speed_range_m_s": (float(speed.min()), float(speed.max())),
            "pressure_range_pa": (float(pressure.min()), float(pressure.max()))}


class _P2Traces:
    """Integrate signed and absolute normal flux along arbitrary straight lines."""

    def __init__(self, fem, state):
        self.fem = fem
        self.vertices = fem.points[fem.triangles]
        self.lower = self.vertices.min(axis=1)
        self.upper = self.vertices.max(axis=1)
        self.bary_inverse = np.linalg.inv((self.vertices[:, 1:] - self.vertices[:, :1]).transpose(0, 2, 1))
        self.velocity = state.values[:2 * fem.n_u].reshape((2, fem.n_u)).T

    def line(self, a, b):
        a, b = np.asarray((a, b), dtype=float) / REFERENCE_HEIGHT_MM
        direction = b - a
        length = np.linalg.norm(direction)
        if length == 0.:
            return 0., 0., 0.
        candidate = np.flatnonzero(np.all(self.upper >= np.minimum(a, b)-1e-12, axis=1)
                                   & np.all(self.lower <= np.maximum(a, b)+1e-12, axis=1))
        inverse = self.bary_inverse[candidate]
        uv0 = np.einsum("eij,ej->ei", inverse, a-self.vertices[candidate, 0])
        delta_uv = inverse @ direction
        alpha = np.column_stack((1.0-uv0.sum(axis=1), uv0))
        beta = np.column_stack((-delta_uv.sum(axis=1), delta_uv))
        positive = beta > 1e-13
        negative = beta < -1e-13
        intersects = np.all(positive | negative | (alpha >= -1e-12), axis=1)
        roots = np.divide(-alpha, beta, out=np.zeros_like(beta), where=positive | negative)
        enter = np.maximum(0., np.max(np.where(positive, roots, -np.inf), axis=1))
        leave = np.minimum(1., np.min(np.where(negative, roots, np.inf), axis=1))
        valid = intersects & ((leave-enter) > 1e-12)
        intervals = np.column_stack((enter[valid], leave[valid]))
        candidate, alpha, beta = candidate[valid], alpha[valid], beta[valid]
        cuts = np.unique(np.round(np.r_[0., intervals.ravel(), 1.], 13))
        forward = backward = covered = 0.
        normal = np.array([direction[1], -direction[0]])
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            if hi - lo <= 1e-12:
                continue
            mid = (lo + hi) / 2
            matches = np.flatnonzero((intervals[:, 0] <= mid) & (mid <= intervals[:, 1]))
            if not len(matches):
                continue
            # On a shared mesh edge either element has the same P2 trace.
            row = matches[0]
            params = np.array([lo, mid, hi])
            barycentric = alpha[row] + params[:, None] * beta[row]
            basis, _ = quadratic_basis(barycentric)
            values = (basis @ self.velocity[self.fem.velocity_dofs[candidate[row]]]) @ normal
            crossing, returning = directional_integrals(*values)
            forward += (hi-lo) * crossing
            backward += (hi-lo) * returning
            covered += hi-lo
        if not np.isclose(covered, 1., rtol=0., atol=1e-8):
            raise ValueError("A nominal interface segment was not covered by the immersed FEM mesh")
        return forward, backward, length * REFERENCE_HEIGHT_MM


def _line_segments(geometry):
    if geometry.is_empty:
        return
    if geometry.geom_type in ("LineString", "LinearRing"):
        coords = list(geometry.coords)
        yield from zip(coords[:-1], coords[1:])
    elif hasattr(geometry, "geoms"):
        for part in geometry.geoms:
            yield from _line_segments(part)


def _flux_stats(values, inlet):
    positive, negative, length = values
    scale = HEIGHT_M * VELOCITY_SCALE
    return {"length_mm": float(length), "positive_m2_per_s": float(positive * scale),
            "negative_m2_per_s": float(negative * scale), "positive_over_inlet": float(positive / inlet),
            "negative_over_inlet": float(negative / inlet), "net_over_inlet": float((positive-negative) / inlet),
            "absolute_over_inlet": float((positive+negative) / inlet)}


def _nominal_solid_leakage(domain, state, grid, geometry, height, inlet):
    from shapely import contains_xy
    from shapely.geometry import LineString, box
    from shapely.geometry.polygon import orient
    from shapely.ops import unary_union

    sharp = box(0., 0., 20., 4.).union(box(9.5, 4., 10.5, 7.))
    envelope = box(grid.x_edges_mm[0], grid.y_edges_mm[0], grid.x_edges_mm[-1], grid.y_edges_mm[-1])
    outer = sharp.union(envelope)
    baffle = box(13., 0., 14., height) if height > 0 else None
    traces = _P2Traces(domain.fem, state)
    outer_walls = outer.boundary
    total = np.zeros(3)
    baffle_faces = np.zeros(3)
    component_boundaries = []
    for index, poly in enumerate(geometry.polygons):
        crossing = np.zeros(3)
        baffle_crossing = np.zeros(3)
        poly = orient(poly, sign=1.)
        for ring in (poly.exterior, *poly.interiors):
            for a, b in _line_segments(ring):
                internal = LineString((a, b)).difference(outer_walls)
                direction = np.subtract(b, a)
                for start, stop in _line_segments(internal):
                    if np.dot(np.subtract(stop, start), direction) < 0.:
                        start, stop = stop, start
                    # Right normal of the oriented fluid polygon points into solid.
                    crossing += traces.line(start, stop)
                if baffle is not None:
                    for start, stop in _line_segments(internal.intersection(baffle.boundary)):
                        if np.dot(np.subtract(stop, start), direction) < 0.:
                            start, stop = stop, start
                        baffle_crossing += traces.line(start, stop)
        total += crossing
        baffle_faces += baffle_crossing
        component_boundaries.append({"fluid_component_index": index,
                                     "internal_boundary": _flux_stats(crossing, inlet),
                                     "baffle_faces": _flux_stats(baffle_crossing, inlet)})
    cut = traces.line((13.5, 0.), (13.5, height)) if baffle is not None else (0., 0., 0.)
    fem = domain.fem
    x, y = (fem.quad_points * REFERENCE_HEIGHT_MM).transpose(2, 0, 1)
    velocity = np.einsum("qa,eai->eqi", fem.basis, traces.velocity[fem.velocity_dofs])
    quadrature_speed = np.linalg.norm(velocity, axis=2)

    def solid_speed_stats(solid):
        # The hard contour can cut a FEM element. Existing volume quadrature
        # samples its indicator, so report sampled and exact solid areas.
        inside = contains_xy(solid, x, y)
        weights = fem.weights[inside]
        speed = quadrature_speed[inside]
        area = weights.sum()
        return {"exact_area_mm2": float(solid.area),
                "sampled_area_mm2": float(area * REFERENCE_HEIGHT_MM**2),
                "sampled_rms_speed_m_s": float(VELOCITY_SCALE*np.sqrt(weights @ speed**2/area)) if area else None,
                "sampled_max_speed_m_s": float(VELOCITY_SCALE*speed.max()) if area else None}

    stats = {
        "method": "Exact signed integration of piecewise P2 traces on nominal straight boundary segments. "
                  "Solid speed metrics use FEM quadrature sampling of the hard indicator.",
        "all_internal_fluid_to_solid": _flux_stats(total, inlet),
        "baffle_fluid_to_solid_faces": _flux_stats(baffle_faces, inlet),
        "fluid_component_boundaries": component_boundaries,
        "baffle_mid_chord_positive_x_cut": _flux_stats(cut, inlet),
        "all_solid_speed": solid_speed_stats(outer.difference(geometry.geometry)),
    }
    if baffle is not None:
        stats["baffle_only_speed"] = solid_speed_stats(baffle)
    return stats


def _write_figures(output, grid, projected, geometry, immersed, fitted):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path as PlotPath
    import matplotlib.tri as mtri
    from shapely.geometry import box
    from shapely.geometry.polygon import orient

    def draw_polygons(ax, polygons, color):
        for polygon in polygons:
            polygon = orient(polygon, sign=1.0)
            coordinates, codes = [], []
            for ring in (polygon.exterior, *polygon.interiors):
                vertices = list(ring.coords)
                coordinates.extend(vertices)
                codes.extend([PlotPath.MOVETO, *([PlotPath.LINETO] * (len(vertices)-2)), PlotPath.CLOSEPOLY])
            ax.add_patch(PathPatch(PlotPath(coordinates, codes), facecolor=color, edgecolor="none"))

    def save(name, render):
        path = output / f"{name}.png"
        fig = plt.figure(figsize=(12, 4.8), dpi=150, facecolor="white")
        try:
            ax = fig.add_axes((0, 0, 1, 1))
            ax.set(xlim=(-0.5, 20.5), ylim=(-0.7, 7.7), aspect="equal")
            ax.set_axis_off()
            render(ax)
            with path.open("xb") as stream:
                fig.savefig(stream, format="png", dpi=150, facecolor="white")
        finally:
            plt.close(fig)
        return path.name

    def draw_density(ax):
        color = plt.get_cmap("viridis")(1.0)
        draw_polygons(ax, (box(0, 0, 6, 4), box(16, 0, 20, 4), box(9.5, 6, 10.5, 7)), color)
        ax.pcolormesh(grid.x_edges_mm, grid.y_edges_mm, projected.reshape(grid.shape),
                      cmap="viridis", vmin=0, vmax=1, shading="flat", antialiased=False)

    def contours(ax, mesh, scalar, cmap, limits):
        lo, hi = limits
        if lo == hi:
            ax.tripcolor(mesh, scalar, cmap=cmap, vmin=lo-1.0, vmax=hi+1.0)
        else:
            ax.tricontourf(mesh, scalar, levels=np.linspace(lo, hi, 49), cmap=cmap)

    initial_geometry = box(0, 0, 20, 4).union(box(9.5, 4, 10.5, 7))
    names = [save("original_sharp_t_junction", lambda ax: draw_polygons(ax, (initial_geometry,), "#17364d")),
             save("hybrid_projected_density", draw_density),
             save("hybrid_reconstructed_fluid", lambda ax: draw_polygons(ax, geometry.polygons, "#17364d"))]
    field_data = []
    for prefix, fields in (("immersed", immersed), ("fitted", fitted)):
        if fields is None:
            continue
        fem, state = fields
        velocity, pressure = _state_fields(fem, state)
        field_data.append((prefix, fem, np.linalg.norm(velocity, axis=1), pressure))
    color_limits = {}
    if field_data:
        for name, index in (("speed_m_s", 2), ("pressure_pa", 3)):
            fields = [entry[index] for entry in field_data]
            color_limits[name] = (min(float(field.min()) for field in fields),
                                  max(float(field.max()) for field in fields))
    for prefix, fem, speed, pressure in field_data:
        a, b, c, ab, bc, ca = fem.velocity_dofs.T
        subdivided = np.vstack([np.column_stack(nodes) for nodes in ((a, ab, ca), (ab, b, bc),
                                                                      (ca, bc, c), (ab, bc, ca))])
        velocity_mesh = mtri.Triangulation(*(fem.velocity_points * REFERENCE_HEIGHT_MM).T, triangles=subdivided)
        pressure_mesh = mtri.Triangulation(*(fem.points * REFERENCE_HEIGHT_MM).T, triangles=fem.triangles)
        names.append(save(f"{prefix}_velocity_magnitude", lambda ax: contours(
            ax, velocity_mesh, speed, "magma", color_limits["speed_m_s"])))
        names.append(save(f"{prefix}_pressure", lambda ax: contours(
            ax, pressure_mesh, pressure, "cividis", color_limits["pressure_pa"])))
    return names, color_limits


def export_endpoint(output, *, grid, projected, baffle_height_mm, domain=None, state=None,
                    fitted_spacings_mm=(), write_pngs=True):
    """Archive topology, hard baffle geometry, physical checks and field evidence."""
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run hybrid reconstruction and fitted diagnostics on Sol")
    if (domain is None) != (state is None):
        raise ValueError("Supply both the immersed domain and state when exporting immersed fields")
    spacings = tuple(float(spacing) for spacing in fitted_spacings_mm)
    if any(not np.isfinite(h) or h <= 0. for h in spacings) or any(a <= b for a, b in zip(spacings[:-1], spacings[1:])):
        raise ValueError("Fitted spacings must be finite, positive and strictly decreasing")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    projected = np.asarray(projected, dtype=float)
    pre, post = _reconstruction_pair(grid, projected, baffle_height_mm)
    before = write_reconstruction(pre, output / "gamma_0p5_without_baffle.json")
    after = write_reconstruction(post, output / "gamma_0p5_minus_hard_baffle.json")
    with (output / "projected_topology_density.npz").open("xb") as stream:
        np.savez_compressed(stream, x_edges_mm=grid.x_edges_mm, y_edges_mm=grid.y_edges_mm,
                            projected_density=projected.reshape(grid.shape))
    fitted_report, fitted_fields = _fitted_checks(post, spacings)
    immersed_fields = None if domain is None else (domain.fem, state)
    report = {
        "sol_job_id": os.environ["SLURM_JOB_ID"], "baffle_extent_mm": [[13., 14.], [0., baffle_height_mm]],
        "exact_reserved_fluid_strip_mm": [[13., 14.], [0., 4.]],
        "nominal_gap_mm": 4.0 - baffle_height_mm, "baffle_area_removed_from_fluid_mm2": pre.area_mm2 - post.area_mm2,
        "before_baffle_area_mm2": pre.area_mm2, "after_baffle_area_mm2": post.area_mm2,
        "geometric_area_feasible": post.area_mm2 <= MAX_AREA_MM2,
        "connected_before_baffle": before["connected_ports"], "connected_after_baffle": after["connected_ports"],
        "post_subtraction_port_components": after["port_components"], "fitted": fitted_report,
        "fitted_passed": bool(fitted_report.get("fitted_diagnostic_passed", False)),
        "projected_topology_density_archive": "projected_topology_density.npz",
    }
    if immersed_fields is not None:
        flux, mass_error, power, split, _ = response_values(domain, state.values)
        report["immersed"] = {**_write_state(output, "immersed_fields", *immersed_fields),
                              "side_share": split, "dimensionless_flux": flux.tolist(), "power_W_per_m": power,
                              "mass_imbalance": mass_error, "nonlinear_residual": state.residual_history[-1]}
        try:
            report["immersed"]["nominal_solid_leakage"] = _nominal_solid_leakage(
                domain, state, grid, post, baffle_height_mm, flux[0])
        except Exception as error:
            report["immersed"]["nominal_solid_leakage_error"] = f"{type(error).__name__}: {error}"
    if fitted_fields is not None:
        report["fitted_fields"] = _write_state(output, "fitted_fields", *fitted_fields)
        report["fitted_fields"]["spacing_mm"] = fitted_report["levels"][-1]["spacing_mm"]
    if write_pngs:
        try:
            names, color_limits = _write_figures(output, grid, projected, post, immersed_fields, fitted_fields)
            report["images"] = names
            report["shared_field_color_limits"] = color_limits
        except Exception as error:
            report["image_error"] = f"{type(error).__name__}: {error}"
        report["image_extent_mm"] = [-0.5, 20.5, -0.7, 7.7]
        report["image_size_pixels"] = [1800, 720]
    with (output / "hybrid_endpoint_report.json").open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return report
