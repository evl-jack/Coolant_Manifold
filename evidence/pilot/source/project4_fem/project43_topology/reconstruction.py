"""Extract and replay the 0.5 contour as an impermeable planar fluid domain."""

from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from ..fem import TaylorHood
from ..junction import boundary_data
from ..project43_power import port_power
from .domain import REFERENCE_HEIGHT_MM


_PORTS_MM = (((0., 0.), (0., 4.)), ((20., 0.), (20., 4.)), ((9.5, 7.), (10.5, 7.)))


@dataclass(frozen=True)
class Reconstruction:
    geometry: object
    polygons: tuple
    port_components: tuple
    contour_segments_mm: tuple

    @property
    def port_component(self):
        ids = set(self.port_components)
        if len(ids) == 1 and None not in ids:
            return self.polygons[ids.pop()]
        return None

    @property
    def area_mm2(self):
        return float(self.geometry.area)


def _crossing(a, b, fa, fb):
    # Canonical endpoint order preserves exact shared coordinates across cells.
    if b < a:
        a, b, fa, fb = b, a, fb, fa
    if fa == 0.:
        return a
    if fb == 0.:
        return b
    fraction = fa / (fa - fb)
    return tuple(x + fraction * (y - x) for x, y in zip(a, b))


def _cell_regions(corners, values):
    """Return polygon rings and contour segments within one bilinear cell.

    Corners are ordered counterclockwise from bottom left. Marching squares
    approximates each continuous bilinear contour branch by a straight chord.
    """
    fluid = values >= 0.
    if fluid.all():
        return [corners], []
    if not fluid.any():
        return [], []

    edges = {}
    for i in range(4):
        j = (i + 1) % 4
        if fluid[i] != fluid[j]:
            edges[i] = _crossing(corners[i], corners[j], values[i], values[j])

    if len(edges) == 2:
        ring = []
        for i, corner in enumerate(corners):
            if fluid[i]:
                ring.append(corner)
            if i in edges:
                ring.append(edges[i])
        return [ring], [tuple(edges.values())]

    # Test the sign at the true bilinear saddle, without a center average.
    f0, f1, f2, f3 = values
    determinant = f0 * f2 - f1 * f3
    coefficient = (f0 - f1) + (f2 - f3)
    if determinant == 0.:
        u = (f0 - f3) / coefficient
        v = (f0 - f1) / coefficient
        x0, y0 = corners[0]
        x1, y1 = corners[2]
        saddle = (x0 + (x1 - x0) * u, y0 + (y1 - y0) * v)
        rings = [[corners[i], edges[i], saddle, edges[(i - 1) % 4]]
                 for i in range(4) if fluid[i]]
        return rings, [(point, saddle) for point in edges.values()]

    connected = (determinant > 0.) == (coefficient > 0.)
    if connected:
        ring = []
        for i, corner in enumerate(corners):
            if fluid[i]:
                ring.append(corner)
            ring.append(edges[i])
        return [ring], [(edges[(i - 1) % 4], edges[i]) for i in range(4) if not fluid[i]]
    rings = [[corners[i], edges[i], edges[(i - 1) % 4]] for i in range(4) if fluid[i]]
    return rings, [(edges[(i - 1) % 4], edges[i]) for i in range(4) if fluid[i]]


def reconstruct(grid, projected):
    """Union the polygonal threshold contour with the three exact fixed-fluid strips."""
    try:
        from shapely.geometry import LineString, Polygon, box
        from shapely.ops import unary_union
    except ImportError as error:
        raise RuntimeError("Fitted reconstruction requires Shapely on Sol") from error

    projected = np.asarray(projected, dtype=float)
    if projected.shape != (grid.size,) or not np.isfinite(projected).all():
        raise ValueError("Expected one finite projected density value per design cell")
    centers = grid.centers_mm.reshape(*grid.shape, 2)
    x = np.r_[grid.x_edges_mm[0], centers[0, :, 0], grid.x_edges_mm[-1]]
    y = np.r_[grid.y_edges_mm[0], centers[:, 0, 1], grid.y_edges_mm[-1]]
    # Duplicate each boundary row and column to extend the density constantly
    # in the direction normal to the design boundary.
    field = np.pad(projected.reshape(grid.shape) - .5, 1, mode="edge")
    pieces = [box(0., 0., 6., 4.), box(16., 0., 20., 4.), box(9.5, 6., 10.5, 7.)]
    contour = []
    for iy, (y0, y1) in enumerate(zip(y[:-1], y[1:])):
        for ix, (x0, x1) in enumerate(zip(x[:-1], x[1:])):
            corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
            values = field[[iy, iy, iy + 1, iy + 1], [ix, ix + 1, ix + 1, ix]]
            rings, segments = _cell_regions(corners, values)
            for ring in rings:
                polygon = Polygon(ring)
                if polygon.area > 0.:
                    if not polygon.is_valid:
                        raise ValueError(f"Invalid threshold region at design patch {(iy, ix)}")
                    pieces.append(polygon)
            contour.extend(segment for segment in segments if segment[0] != segment[1])

    geometry = unary_union(pieces)
    if not geometry.is_valid or geometry.geom_type not in ("Polygon", "MultiPolygon"):
        raise ValueError("Fixed-fluid union did not produce a valid polygonal geometry")
    polygons = (geometry,) if geometry.geom_type == "Polygon" else tuple(geometry.geoms)
    polygons = tuple(sorted(polygons, key=lambda polygon: (*polygon.bounds, polygon.area)))
    port_components = []
    for port in _PORTS_MM:
        line = LineString(port)
        matches = [i for i, polygon in enumerate(polygons) if polygon.boundary.covers(line)]
        if len(matches) != 1:
            raise ValueError(f"The exact physical port {port} was not preserved")
        port_components.append(matches[0])
    return Reconstruction(geometry, polygons, tuple(port_components), tuple(contour))


def write_reconstruction(result, path):
    """Archive threshold segments, final polygons, islands and port connectivity."""
    from shapely.geometry import mapping

    geometry = mapping(result.geometry)
    digest = sha256(json.dumps(geometry, sort_keys=True).encode()).hexdigest()
    data = dict(units="mm", threshold=.5, geometry=geometry, geometry_sha256=digest,
                raw_threshold_segments_mm=result.contour_segments_mm,
                physical_fluid_area_mm2=result.area_mm2,
                port_order=("inlet", "main", "side"), port_components=result.port_components,
                connected_ports=result.port_component is not None,
                components=[dict(area_mm2=p.area, solid_islands=len(p.interiors), bounds_mm=p.bounds)
                            for p in result.polygons])
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return data


def fitted_mesh(polygon, max_spacing_mm):
    """Use a boundary-edge spacing and interior triangle-area target in mm."""
    try:
        import triangle
    except ImportError as error:
        raise RuntimeError("Fitted mesh generation requires Triangle on Sol") from error
    from shapely.geometry import Polygon

    if not np.isfinite(max_spacing_mm) or max_spacing_mm <= 0.:
        raise ValueError("Fitted spacing must be finite and positive")
    points = []
    point_ids = {}
    segments = []
    holes = []
    for ring in (polygon.exterior, *polygon.interiors):
        original = np.asarray(ring.coords[:-1])
        loop = []
        for a, b in zip(original, np.roll(original, -1, axis=0)):
            count = max(1, int(np.ceil(np.linalg.norm(b - a) / max_spacing_mm)))
            loop.extend(np.linspace(a, b, count, endpoint=False))
        ids = []
        for point in loop:
            key = tuple(point)
            if key not in point_ids:
                point_ids[key] = len(points)
                points.append(point)
            ids.append(point_ids[key])
        segments.extend((a, b) for a, b in zip(ids, np.roll(ids, -1)) if a != b)
    for ring in polygon.interiors:
        holes.append(tuple(Polygon(ring).representative_point().coords[0]))

    pslg = {"vertices": np.asarray(points), "segments": np.asarray(segments)}
    if holes:
        pslg["holes"] = np.asarray(holes)
    max_area = .5 * max_spacing_mm ** 2
    mesh = triangle.triangulate(pslg, f"pq25a{max_area:.16f}Q")
    if "triangles" not in mesh:
        raise ValueError("Triangle returned no fitted elements")
    points = mesh["vertices"]
    tris = mesh["triangles"].copy()
    ab = points[tris[:, 1]] - points[tris[:, 0]]
    ac = points[tris[:, 2]] - points[tris[:, 0]]
    twice_area = ab[:, 0] * ac[:, 1] - ab[:, 1] * ac[:, 0]
    flipped = twice_area < 0.
    tris[flipped] = tris[flipped][:, [0, 2, 1]]
    np.testing.assert_allclose(np.abs(twice_area).sum() / 2, polygon.area, rtol=1e-10)
    ids, compact = np.unique(tris, return_inverse=True)
    return points[ids] / REFERENCE_HEIGHT_MM, compact.reshape(-1, 3)


def solve_fitted(result, max_spacing_mm):
    """Solve the port-connected polygon and account separately for closed pockets."""
    if result.port_component is None:
        raise ValueError("No reconstructed component connects all three ports")
    fem = TaylorHood(*fitted_mesh(result.port_component, max_spacing_mm))
    fixed, prescribed, ports, flux = boundary_data(fem)
    widths = []
    for edges, expected in zip(ports, _PORTS_MM):
        a, b = np.asarray(expected)
        endpoints = fem.points[edges[:, :2]] * REFERENCE_HEIGHT_MM
        np.testing.assert_allclose(endpoints.min(axis=(0, 1)), a, rtol=0., atol=1e-12)
        np.testing.assert_allclose(endpoints.max(axis=(0, 1)), b, rtol=0., atol=1e-12)
        width = np.linalg.norm(endpoints[:, 1] - endpoints[:, 0], axis=1).sum()
        np.testing.assert_allclose(width, np.linalg.norm(b - a), rtol=1e-12)
        widths.append(width)

    state = None
    for reynolds in (0., 25., 50.):
        fem.inertia = reynolds
        state = fem.solve(fixed, prescribed, fem.load(), initial=None if state is None else state.values)
    inlet, main, side = flux @ state.values
    mass_error = abs(inlet - main - side) / inlet
    speed = 50. * .001002 / (998.2 * .004)
    power = sum(port_power(fem, ports, state.values, 50.)[0]) * .001002 * speed ** 2
    report = dict(triangles=len(fem.triangles), dofs=fem.size, spacing_mm=max_spacing_mm,
                  physical_fluid_area_mm2=result.area_mm2,
                  port_component_area_mm2=result.port_component.area,
                  closed_pocket_area_mm2=result.area_mm2 - result.port_component.area,
                  closed_pocket_count=len(result.polygons) - 1, port_widths_mm=widths,
                  side_share=side / inlet, inlet_flux=inlet, main_flux=main, side_flux=side,
                  mechanical_power_W_per_m=power, mass_imbalance=mass_error,
                  nonlinear_residual=state.residual_history[-1],
                  state_passed=bool(min(inlet, main, side) > 0. and mass_error <= 1e-9
                                    and state.residual_history[-1] <= 1e-10))
    return report, fem, state
