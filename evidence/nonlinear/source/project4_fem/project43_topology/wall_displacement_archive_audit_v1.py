"""Read completed Sol wall-displacement checkpoints and independently audit their geometry."""

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path

import numpy as np
from shapely.geometry import LineString, Point, Polygon, box, shape
from shapely.ops import unary_union

from .density import DesignGrid
from .reconstruction import reconstruct


PORTS = (((0., 0.), (0., 4.)), ((20., 0.), (20., 4.)), ((9.5, 7.), (10.5, 7.)))
CORNERS = ((9.5, 4.), (10.5, 4.))
SECANT_MM, SAMPLE_MM, LIMIT_DEG = .25, .025, 45.
TOL = 1e-10


def _read(path):
    contents = path.read_bytes()
    doc = json.loads(contents)
    return doc.get("results", doc), sha256(contents).hexdigest()


def _cell_footprints(grid):
    x, y = grid.x_edges_mm, grid.y_edges_mm
    return [box(a, c, b, d) for c, d in zip(y[:-1], y[1:]) for a, b in zip(x[:-1], x[1:])]


def _independent_mask(cells):
    window = box(6.5, 0., 15.5, 5.5)
    core = box(0., 1., 20., 3.).union(box(9.5, 2., 10.5, 7.))
    walls = tuple(LineString(pair) for pair in (
        ((0., 0.), (20., 0.)), ((0., 4.), (9.5, 4.)), ((10.5, 4.), (20., 4.)),
        ((9.5, 4.), (9.5, 7.)), ((10.5, 4.), (10.5, 7.)),
    ))
    return np.array([window.covers(cell) and cell.intersection(core).area == 0.
                     and any(all(wall.distance(Point(p)) <= 1. for p in cell.exterior.coords[:-1])
                             for wall in walls) for cell in cells])


def _anchors_and_filter(grid, free, raw, filtered, projected, design, beta):
    fields = (raw, filtered, projected)
    for field in fields:
        if field.shape != (grid.size,) or not np.isfinite(field).all():
            raise ValueError("A full archived density field has an unexpected shape or nonfinite values")
    f = filtered.reshape(grid.shape)
    lap = np.zeros_like(f)
    for axis in (0, 1):
        left = [slice(None), slice(None)]
        right = left.copy()
        left[axis], right[axis] = slice(None, -1), slice(1, None)
        left, right = tuple(left), tuple(right)
        delta = f[left] - f[right]
        lap[left] += delta
        lap[right] -= delta
    residual = grid.cell_area_mm2 * (filtered - raw) + .5 ** 2 * lap.ravel()
    replay = filtered if beta == 0. else .5 + np.tanh(beta * (filtered - .5)) / (2. * np.tanh(beta / 2.))
    return dict(compact_design_equals_raw_editable=bool(np.array_equal(design, raw[free])),
                inactive_anchors_exact={name: bool(np.array_equal(field[~free], grid.base_raw[~free]))
                                        for name, field in zip(("raw", "filtered", "projected"), fields)},
                filter_residual_inf=float(np.max(np.abs(residual[free]))),
                projection_abs_error_inf=float(np.max(np.abs(replay - projected))))


def _ring_data(ring):
    vertices = np.array(ring.coords)
    lengths = np.r_[0., np.cumsum(np.linalg.norm(np.diff(vertices, axis=0), axis=1))]
    return LineString(vertices), vertices, lengths


def _angle(vertices, lengths, s):
    arcs = (np.asarray([s - SECANT_MM, s, s + SECANT_MM])) % lengths[-1]
    a, b, c = np.column_stack([np.interp(arcs, lengths, vertices[:, i]) for i in range(2)])
    u, v = b - a, c - b
    if lengths[-1] <= 2 * SECANT_MM or np.linalg.norm(u) * np.linalg.norm(v) <= TOL ** 2:
        return b, None
    determinant = u[0] * v[1] - u[1] * v[0]
    return b, float(np.degrees(np.arctan2(abs(determinant), u @ v)))


def _turns(polygons, footprint):
    ports = [LineString(ends) for ends in PORTS]
    rings = [(i, j, *_ring_data(ring)) for i, p in enumerate(polygons)
             for j, ring in enumerate((p.exterior, *p.interiors))]
    turns = []

    def add(ring, s):
        component, index, _, vertices, lengths = ring
        point, angle = _angle(vertices, lengths, s)
        sample = dict(component=component, ring=index, arc_mm=float(s), xy_mm=point.tolist(), angle_degrees=angle)
        turns.append(sample)
        return sample

    for ring in rings:
        _, _, line, _, lengths = ring
        for s in np.unique(np.r_[lengths[:-1], np.arange(0., lengths[-1], SAMPLE_MM)]):
            point = line.interpolate(float(s))
            if point.distance(footprint) <= SECANT_MM + TOL and all(point.distance(port) > TOL for port in ports):
                add(ring, s)
    junctions = []
    for xy in CORNERS:
        point = Point(xy)
        ring = min(rings, key=lambda r: r[2].distance(point))
        junctions.append(dict(original_xy_mm=xy, wall_distance_mm=float(ring[2].distance(point)),
                              nearest_wall_sample=add(ring, ring[2].project(point))))
    finite = [s for s in turns if s["angle_degrees"] is not None]
    maximum = max(finite, key=lambda s: s["angle_degrees"]) if finite else None
    undefined = len(turns) - len(finite)
    failures = [s for s in finite if s["angle_degrees"] > LIMIT_DEG + TOL]
    return dict(turning_angle_passed=bool(maximum and not failures and not undefined), max_turn=maximum,
                sample_count=len(turns), over_limit_count=len(failures), undefined_count=undefined,
                original_junction_corner_samples=junctions)


def _audit_checkpoint(folder, grid, cells, free, original):
    archived, checkpoint_sha = _read(folder / "checkpoint.json")
    report = dict(checkpoint_sha256=checkpoint_sha, stage=archived.get("stage"),
                  qualified_optimization_terminal=archived.get("qualified_optimization_terminal"),
                  side_fraction=archived.get("side_fraction"), power_w_per_m=archived.get("power_w_per_m"))
    if archived.get("error"):
        return dict(report, archived_error=archived["error"])
    fields_path = folder / "fields.npz"
    keys = (("raw", "raw_density"), ("filtered", "filtered_density"),
            ("projected", "projected_density"), ("compact_design", "design"))
    with np.load(fields_path, allow_pickle=False) as fields:
        raw, filtered, projected, design = (fields[key] for key, _ in keys)
        matches = {key: bool(np.array_equal(fields[key], archived[name])) for key, name in keys}
        state = fields["state"]
        state_finite = bool(state.ndim == 1 and state.size and np.isfinite(state).all())
    report.update(fields_sha256=sha256(fields_path.read_bytes()).hexdigest(), fields_match_checkpoint=matches,
                  state_finite=state_finite, **_anchors_and_filter(
                      grid, free, raw, filtered, projected, design, float(archived["stage"]["beta"])))
    if not (folder / "geometry.json").is_file():
        return dict(report, archived_geometry_missing=True, archived_geometry_audit=archived.get("wall_geometry_audit"))
    geometry, geometry_sha = _read(folder / "geometry.json")
    final = shape(geometry["geometry"])
    if final.is_empty or not final.is_valid or final.geom_type not in ("Polygon", "MultiPolygon"):
        raise ValueError("Archived fluid geometry is not a valid polygon or multipolygon")
    polygons = [final] if final.geom_type == "Polygon" else list(final.geoms)
    cells = [cell for cell, editable in zip(cells, free) if editable]
    footprint = unary_union(cells)
    raw_shape = reconstruct(grid, projected).geometry
    inside = unary_union([raw_shape.intersection(cell) for cell in cells])
    expected = original.difference(footprint).union(inside)
    outside_error = original.symmetric_difference(final).difference(footprint).area
    formula_error = expected.symmetric_difference(final).area
    port_indices = [tuple(i for i, p in enumerate(polygons) if p.boundary.covers(LineString(pair))) for pair in PORTS]
    pocket_areas = [float(p.area) for i, p in enumerate(polygons) if not any(i in ids for ids in port_indices)]
    hole_areas = [float(Polygon(ring).area) for p in polygons for ring in p.interiors]
    turns = _turns(polygons, footprint)
    archived_turn = (archived.get("wall_geometry_audit") or {}).get("max_turn_degrees")
    return dict(report, geometry_file_sha256=geometry_sha, geometry_matches_checkpoint=geometry == archived.get("geometry"),
                shape_digest_matches=sha256(json.dumps(geometry["geometry"], sort_keys=True).encode()).hexdigest()
                == geometry.get("geometry_sha256"), exact_port_segments_preserved=all(len(ids) == 1 for ids in port_indices),
                connected_ports=len(port_indices[0]) == 1 and port_indices[0] == port_indices[1] == port_indices[2],
                closed_pocket_count=len(pocket_areas), closed_pocket_areas_mm2=pocket_areas,
                hole_count=len(hole_areas), hole_areas_mm2=hole_areas, fluid_area_mm2=float(final.area),
                outside_editable_difference_mm2=float(outside_error), clip_formula_difference_mm2=float(formula_error),
                exact_clipping_passed=max(outside_error, formula_error) <= TOL,
                archived_max_turn_degrees=archived_turn, **turns)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    job = os.environ.get("SLURM_JOB_ID")
    if not job:
        parser.error("Run numerical archive audits on Sol under Slurm")
    setup, setup_sha = _read(args.run / "setup.json")
    grid = DesignGrid(setup["design_spacing_mm"])
    cells = _cell_footprints(grid)
    free = _independent_mask(cells)
    if not np.array_equal(free, setup["free_cells"]):
        raise ValueError("Archived editable cells do not match the independently derived mask")
    original = box(0., 0., 20., 4.).union(box(9.5, 4., 10.5, 7.))
    results = []
    files = sorted(args.run.glob("stage_*/*/checkpoint.json"))
    for path in files:
        row = {"checkpoint": path.relative_to(args.run).as_posix()}
        try:
            row.update(_audit_checkpoint(path.parent, grid, cells, free, original))
        except (OSError, ValueError, TypeError, KeyError, RuntimeError) as error:
            row["read_or_audit_error"] = f"{type(error).__name__}: {error}"
        results.append(row)
    result = dict(input_dir=str(args.run), archive_job_id=setup.get("slurm_job_id"), audit_job_id=job,
                  setup_sha256=setup_sha, independently_verified_editable_cell_count=int(free.sum()),
                  secant_length_mm=SECANT_MM, arc_step_mm=SAMPLE_MM, turn_limit_degrees=LIMIT_DEG,
                  completed_checkpoint_files=len(files), results=results)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
