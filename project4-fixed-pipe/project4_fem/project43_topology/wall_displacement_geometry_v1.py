"""Exact edit-footprint clipping and wall-turn screening for wall displacement."""

import numpy as np
from shapely.geometry import LineString, Point, Polygon, box
from shapely.ops import unary_union

from .reconstruction import Reconstruction, _PORTS_MM, reconstruct


JUNCTION_CORNERS_MM = ((9.5, 4.), (10.5, 4.))
SECANT_LENGTH_MM = .25
ARC_STEP_MM = .025
TURN_LIMIT_DEGREES = 45.
AREA_ROUNDOFF_MM2 = 1e-10
COORDINATE_ROUNDOFF_MM = 1e-10


def sharp_geometry():
    return box(0., 0., 20., 4.).union(box(9.5, 4., 10.5, 7.))


def editable_geometry(grid, free):
    """Union the exact square footprints of the editable, row-major cells."""
    free = np.asarray(free)
    if free.dtype.kind == "b":
        if free.shape not in ((grid.size,), grid.shape):
            raise ValueError("The editable mask must cover the entire design grid")
        free = np.flatnonzero(free)
    if (free.ndim != 1 or free.dtype.kind not in "iu" or free.size == 0
            or np.any(free < 0) or np.any(free >= grid.size) or len(np.unique(free)) != free.size):
        raise ValueError("Expected distinct, valid editable design-cell IDs")
    iy, ix = np.unravel_index(free.astype(np.intp, copy=False), grid.shape)
    x, y = grid.x_edges_mm, grid.y_edges_mm
    return unary_union([box(x[i], y[j], x[i + 1], y[j + 1]) for j, i in zip(iy, ix)])


def _polygon_parts(shape):
    if shape.is_empty:
        return []
    if shape.geom_type == "Polygon":
        return [shape]
    if shape.geom_type in ("MultiPolygon", "GeometryCollection"):
        return [polygon for part in shape.geoms for polygon in _polygon_parts(part)]
    return []


def _clipped(original, editable, raw):
    final = unary_union(_polygon_parts(original.difference(editable))
                        + _polygon_parts(raw.geometry.intersection(editable)))
    if final.is_empty or not final.is_valid or final.geom_type not in ("Polygon", "MultiPolygon"):
        raise ValueError("Editable-cell clipping did not produce valid polygonal fluid geometry")
    parts = tuple(sorted(_polygon_parts(final), key=lambda p: (*p.bounds, p.area)))
    ports = []
    for endpoints in _PORTS_MM:
        port = LineString(endpoints)
        matches = [i for i, polygon in enumerate(parts) if polygon.boundary.covers(port)]
        if len(matches) != 1:
            raise ValueError(f"Clipping changed exact physical port {endpoints}")
        ports.append(matches[0])
    return Reconstruction(final, parts, tuple(ports), raw.contour_segments_mm)


def reconstruct_clipped(grid, projected, free):
    """Keep the sharp geometry outside editable cells and threshold density only inside."""
    return _clipped(sharp_geometry(), editable_geometry(grid, free), reconstruct(grid, projected))


def _wall_samples(polygon):
    for ring_id, ring in enumerate((polygon.exterior, *polygon.interiors)):
        wall = LineString(ring.coords)
        points = np.asarray(wall.coords)
        vertices = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))][:-1]
        steps = np.arange(0., wall.length, ARC_STEP_MM)
        positions = np.unique(np.r_[vertices, steps])
        yield ring_id, wall, positions


def _turn(wall, position):
    arc = (np.asarray([position - SECANT_LENGTH_MM, position, position + SECANT_LENGTH_MM])
           % wall.length)
    a, b, c = (np.asarray(wall.interpolate(float(s)).coords[0]) for s in arc)
    incoming, outgoing = b - a, c - b
    count = np.linalg.norm(incoming) * np.linalg.norm(outgoing)
    if wall.length <= 2 * SECANT_LENGTH_MM or count <= COORDINATE_ROUNDOFF_MM ** 2:
        return b, None
    cosine = np.clip(np.dot(incoming, outgoing) / count, -1., 1.)
    return b, float(np.degrees(np.arccos(cosine)))


def _report_sample(ring, arc, xy, degrees):
    return dict(ring=ring, arc_mm=float(arc), xy_mm=tuple(float(x) for x in xy), degrees=degrees)


def _turning_report(original, final, editable):
    # At a clipping join, a secant can cross the edit footprint even when its
    # center lies just outside it. Include every center within one secant span.
    ports = tuple(LineString(port) for port in _PORTS_MM)
    port_tips = tuple(Point(xy) for endpoints in _PORTS_MM for xy in endpoints)
    original_wall = LineString(original.exterior.coords)
    corners = []
    for xy in JUNCTION_CORNERS_MM:
        point = Point(xy)
        arc = original_wall.project(point)
        _, baseline_degrees = _turn(original_wall, arc)
        corners.append(dict(original_xy_mm=xy, original_turn_degrees=baseline_degrees))

    walls = []
    checked = []
    undefined = []
    for component, polygon in enumerate(final.polygons):
        for ring, wall, positions in _wall_samples(polygon):
            walls.append((component, ring, wall))
            for arc in positions:
                point = wall.interpolate(float(arc))
                if point.distance(editable) > SECANT_LENGTH_MM + COORDINATE_ROUNDOFF_MM:
                    continue
                if any(point.distance(port) <= COORDINATE_ROUNDOFF_MM for port in ports):
                    continue
                if any(point.distance(tip) <= COORDINATE_ROUNDOFF_MM for tip in port_tips):
                    continue
                xy, degrees = _turn(wall, arc)
                sample = dict(component=component, **_report_sample(ring, arc, xy, degrees))
                if degrees is None:
                    undefined.append(sample)
                else:
                    checked.append(sample)

    for entry in corners:
        point = Point(entry["original_xy_mm"])
        component, ring, wall = min(walls, key=lambda item: item[2].distance(point))
        arc = wall.project(point)
        xy, degrees = _turn(wall, arc)
        sample = dict(component=component, **_report_sample(ring, arc, xy, degrees))
        entry.update(final_wall_distance_mm=float(wall.distance(point)), final_sample=sample)
        if degrees is None:
            undefined.append(sample)
        else:
            checked.append(sample)

    max_sample = max(checked, key=lambda sample: sample["degrees"]) if checked else None
    maximum = max_sample["degrees"] if max_sample else None
    return dict(secant_length_mm=SECANT_LENGTH_MM, arc_step_mm=ARC_STEP_MM,
                limit_degrees=TURN_LIMIT_DEGREES, max_turn_degrees=maximum,
                max_turn_sample=max_sample, sample_count=len(checked),
                above_limit_count=sum(s["degrees"] > TURN_LIMIT_DEGREES + 1e-10 for s in checked),
                undefined_samples=undefined, original_junction_corners=corners,
                turning_angle_passed=bool(not undefined and maximum is not None
                                          and maximum <= TURN_LIMIT_DEGREES + 1e-10))


def audit_wall_displacement(grid, projected, free, result=None):
    """Report exact clipping, topology, area and sampled wall angles."""
    original = sharp_geometry()
    editable = editable_geometry(grid, free)
    raw = reconstruct(grid, projected)
    final = _clipped(original, editable, raw) if result is None else result
    expected = original.difference(editable).union(raw.geometry.intersection(editable))
    outside_error = original.symmetric_difference(final.geometry).difference(editable).area
    clip_error = final.geometry.symmetric_difference(expected).area
    connected = final.port_component is not None
    pockets = [p for i, p in enumerate(final.polygons) if i not in final.port_components]
    holes = [Polygon(ring) for p in final.polygons for ring in p.interiors]
    report = dict(connected_ports=connected, exact_port_segments_preserved=all(
                      sum(p.boundary.covers(LineString(endpoints)) for p in final.polygons) == 1
                      for endpoints in _PORTS_MM),
                  port_components=final.port_components, component_count=len(final.polygons),
                  closed_pocket_count=len(pockets), closed_pocket_areas_mm2=[float(p.area) for p in pockets],
                  hole_count=len(holes), hole_areas_mm2=[float(p.area) for p in holes],
                  fluid_area_mm2=final.area_mm2, sharp_fluid_area_mm2=float(original.area),
                  main_port_component_area_mm2=float(final.port_component.area) if connected else None,
                  editable_area_mm2=float(editable.area),
                  editable_cell_count=int(round(editable.area / grid.cell_area_mm2)),
                  added_fluid_area_mm2=float(final.geometry.difference(original).area),
                  removed_fluid_area_mm2=float(original.difference(final.geometry).area),
                  raw_fluid_area_mm2=raw.area_mm2,
                  raw_outside_editable_symmetric_difference_mm2=float(
                      original.symmetric_difference(raw.geometry).difference(editable).area),
                  outside_editable_symmetric_difference_mm2=float(outside_error),
                  clip_formula_symmetric_difference_mm2=float(clip_error),
                  exact_clipping_passed=bool(max(outside_error, clip_error) <= AREA_ROUNDOFF_MM2))
    report.update(_turning_report(original, final, editable))
    return report
