"""Export matched original-FEM evidence for archived baffle-free manifold cases."""

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path
import re

import numpy as np

from .density import DesignGrid, DensityPipeline
from .domain import REFERENCE_HEIGHT_MM
from .evaluator import HEIGHT_M, VELOCITY_SCALE, VISCOSITY
from .reconstruction import Reconstruction, reconstruct, solve_fitted, write_reconstruction


VIEW_MM = (-0.5, 20.5, -0.7, 7.7)
DEFAULT_FITTED_SPACING_MM = 0.125
FLOW_SPLIT_SCALE = 0.001


def _read(path):
    wrapped = json.loads(path.read_text(encoding="utf-8"))
    if not wrapped.get("completed"):
        raise ValueError(f"The archived report did not complete: {path}")
    return wrapped["results"]


def _write(path, data):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write("\n")


def _load_case(label, lane, checkpoint):
    setup = _read(lane / "setup.json")
    summary_file = lane / "summary.json"
    summary = _read(summary_file) if summary_file.is_file() else {}
    if checkpoint is None:
        best = summary.get("best_physical")
        selected = summary.get("completed_stage_selections", [])
        if best is not None:
            checkpoint = lane / best["checkpoint"] / "checkpoint.json"
        elif selected:
            checkpoint = lane / selected[-1] / "checkpoint.json"
        else:
            selected_files = sorted(lane.glob("stage_*/selected.json"))
            if selected_files:
                checkpoint = lane / _read(selected_files[-1])["checkpoint"] / "checkpoint.json"
    case = {"label": label, "lane": str(lane), "setup": setup, "archive": None, "checkpoint": None}
    if checkpoint is None:
        case["reason"] = "The lane has no archived physical or immersed candidate"
        return case
    archived = _read(checkpoint)
    case.update(archive=str(checkpoint), checkpoint=archived)
    if "design" not in archived:
        raise ValueError(f"No archived design vector at {checkpoint}")
    stage = archived.get("stage")
    if stage is None:
        stage_dir = next(part for part in checkpoint.parts if re.fullmatch(r"stage_[0-9]+", part))
        stage = setup["stages"][int(stage_dir.split("_")[1])]
    case["stage"] = stage
    return case


def _history(case):
    lane = Path(case["lane"])
    setup = case["setup"]
    reference = float(setup["normalization_W_per_m"])
    lower, upper = setup["flow_band"]
    samples = []
    for path in sorted(lane.glob("stage_*/mma.json")):
        archived = _read(path)
        stage = int(path.parent.name.split("_")[1])
        for row in archived["result"]["evaluation_history"]:
            side_share = lower - FLOW_SPLIT_SCALE * row["inequalities"][0]
            upper_replay = upper + FLOW_SPLIT_SCALE * row["inequalities"][1]
            if not np.isclose(side_share, upper_replay, rtol=0., atol=1e-10):
                raise ValueError(f"Inconsistent split history at {path}")
            samples.append({"evaluation": len(samples) + 1, "stage": stage, "side_flow_pct": 100 * side_share,
                            "power_uW_per_m": 1e6 * reference * row["objective"],
                            "immersed_constraints_passed": bool(row["feasible"])})
    return samples


def _latest_common_spacing(cases):
    fitted_levels = []
    for case in cases:
        if case["checkpoint"] is not None:
            spacings = {float(row["spacing_mm"]) for row in case["checkpoint"].get("fitted", {}).get("levels", [])}
            if spacings:
                fitted_levels.append(spacings)
    return min(set.intersection(*fitted_levels), default=DEFAULT_FITTED_SPACING_MM) if fitted_levels else DEFAULT_FITTED_SPACING_MM


def _state_archive(path, fem, state):
    values = np.asarray(state.values)
    velocity = VELOCITY_SCALE * values[:2 * fem.n_u].reshape((2, fem.n_u)).T
    pressure = VISCOSITY * VELOCITY_SCALE / HEIGHT_M * values[2 * fem.n_u:]
    with path.open("xb") as stream:
        np.savez_compressed(stream, vertices_mm=fem.points * REFERENCE_HEIGHT_MM, triangles=fem.triangles,
                            velocity_vertices_mm=fem.velocity_points * REFERENCE_HEIGHT_MM,
                            velocity_dofs=fem.velocity_dofs, velocity_m_s=velocity, pressure_pa=pressure)
    return fem, np.linalg.norm(velocity, axis=1), pressure


def _fitted_replay(case, geometry, spacing, output):
    if geometry.port_component is None:
        return {"attempted": False, "reason": "Reconstructed ports are disconnected"}, None
    if geometry.area_mm2 > 83.:
        return {"attempted": False, "reason": "Reconstructed fluid area exceeds the cap"}, None
    try:
        metrics, fem, state = solve_fitted(geometry, spacing)
    except Exception as error:
        return {"attempted": True, "error": f"{type(error).__name__}: {error}"}, None
    result = dict(metrics, attempted=True)
    rows = (case.get("checkpoint") or {}).get("fitted", {}).get("levels", [])
    matching = next((row for row in rows if float(row["spacing_mm"]) == spacing), None)
    if matching is not None:
        keys = ("side_share", "mechanical_power_W_per_m", "physical_fluid_area_mm2")
        signed_differences = {name: metrics[name] - matching[name] for name in keys}
        limit = {"side_share": 1e-9, "mechanical_power_W_per_m": 1e-12, "physical_fluid_area_mm2": 1e-8}
        result["matched_archive_differences"] = signed_differences
        result["matched_archive_verified"] = all(abs(value) <= limit[key] for key, value in signed_differences.items())
    return result, _state_archive(output / "fitted_fields.npz", fem, state)


def _candidate_geometry(grid, projected, archived):
    from shapely.geometry import shape

    computed = reconstruct(grid, projected)
    reference = archived.get("geometry", {})
    if not reference.get("geometry_sha256") or not reference.get("geometry"):
        return computed, None
    geometry = reference["geometry"]
    digest = sha256(json.dumps(geometry, sort_keys=True).encode()).hexdigest()
    if digest != reference["geometry_sha256"]:
        raise ValueError("The archived geometry digest did not verify")
    source = shape(geometry)
    area_error = float(computed.geometry.symmetric_difference(source).area)
    if area_error > 1e-10:
        raise ValueError(f"The reconstructed density and archived polygon differ by {area_error:.3e} mm²")
    # Archive coordinates can order polygons differently across GEOS versions.
    # Preserve the polygon coordinates from the endpoint's original report.
    polygons = (source,) if source.geom_type == "Polygon" else tuple(source.geoms)
    polygons = tuple(sorted(polygons, key=lambda polygon: (*polygon.bounds, polygon.area)))
    return Reconstruction(source, polygons, tuple(reference["port_components"]), computed.contour_segments_mm), area_error


def _graphics(root, grid, visuals, histories):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import PathPatch
    from matplotlib.tri import Triangulation
    from shapely.geometry import box
    from shapely.geometry.polygon import orient

    def polygons(ax, shapes, *, color):
        for poly in shapes:
            vertices, codes = [], []
            polygon = orient(poly, sign=1.)
            for ring in (polygon.exterior, *polygon.interiors):
                coords = list(ring.coords)
                vertices.extend(coords)
                codes.extend([MplPath.MOVETO, *([MplPath.LINETO] * (len(coords) - 2)), MplPath.CLOSEPOLY])
            ax.add_patch(PathPatch(MplPath(vertices, codes), facecolor=color, edgecolor="none"))

    def save_frame(folder, name, render):
        path = folder / f"{name}.png"
        fig = plt.figure(figsize=(12, 4.8), dpi=150, facecolor="white")
        try:
            ax = fig.add_axes((0., 0., 1., 1.))
            ax.set(aspect="equal", xlim=VIEW_MM[:2], ylim=VIEW_MM[2:])
            ax.set_axis_off()
            render(ax)
            with path.open("xb") as stream:
                fig.savefig(stream, format="png", dpi=150, facecolor="white")
        finally:
            plt.close(fig)
        return str(path.relative_to(root))

    def render_density(ax, visual):
        collars = (box(0, 0, 6, 4), box(16, 0, 20, 4), box(9.5, 6, 10.5, 7))
        polygons(ax, collars, color=density_cmap(1.))
        ax.pcolormesh(grid.x_edges_mm, grid.y_edges_mm, visual["density"].reshape(grid.shape),
                      cmap=density_cmap, vmin=0, vmax=1, shading="flat", antialiased=False)

    def contours(ax, mesh, values, limits, cmap):
        if limits[0] == limits[1]:
            ax.tripcolor(mesh, values, cmap=cmap, vmin=limits[0] - 1., vmax=limits[1] + 1.)
        else:
            ax.tricontourf(mesh, values, levels=np.linspace(*limits, 73), cmap=cmap)

    fields = [visual["fields"] for visual in visuals if visual["fields"] is not None]
    limits = {"density": (0., 1.)}
    if fields:
        limits["velocity_m_s"] = (0., max(float(field[1].max()) for field in fields))
        limits["pressure_pa"] = (min(float(field[2].min()) for field in fields),
                                  max(float(field[2].max()) for field in fields))
    density_cmap = plt.get_cmap("viridis")
    images = []
    for visual in visuals:
        folder = visual["folder"]
        images.append(save_frame(folder, "geometry", lambda ax: polygons(
            ax, visual["geometry"].polygons, color="#17415c")))
        images.append(save_frame(folder, "density", lambda ax: render_density(ax, visual)))
        if visual["fields"] is None:
            continue
        fem, speed, pressure = visual["fields"]
        a, b, c, ab, bc, ca = fem.velocity_dofs.T
        p2_subtriangles = np.vstack([np.column_stack(nodes) for nodes in (
            (a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca))])
        velocity_mesh = Triangulation(*(fem.velocity_points * REFERENCE_HEIGHT_MM).T, triangles=p2_subtriangles)
        pressure_mesh = Triangulation(*(fem.points * REFERENCE_HEIGHT_MM).T, triangles=fem.triangles)
        images.append(save_frame(folder, "fitted_velocity_magnitude", lambda ax: contours(
            ax, velocity_mesh, speed, limits["velocity_m_s"], "magma")))
        images.append(save_frame(folder, "fitted_pressure", lambda ax: contours(
            ax, pressure_mesh, pressure, limits["pressure_pa"], "cividis")))

    for key, title, filename in (("side_flow_pct", "Side-outlet flow (%)", "flow_history.png"),
                                  ("power_uW_per_m", "Pumping power (µW/m)", "power_history.png")):
        if not any(histories.values()):
            break
        fig, ax = plt.subplots(figsize=(12, 4.8), dpi=150, layout="constrained")
        palette = ("#1b577a", "#b34f25", "#28764d", "#7b5096", "#816827")
        for series, (label, history) in enumerate(histories.items()):
            color = palette[series % len(palette)]
            stages = sorted({row["stage"] for row in history})
            for index, stage in enumerate(stages):
                rows = [row for row in history if row["stage"] == stage]
                x = [row["evaluation"] for row in rows]
                y = [row[key] for row in rows]
                ax.plot(x, y, color=color, linewidth=1.2, alpha=.65, label=label if index == 0 else None)
                feasible = [row for row in rows if row["immersed_constraints_passed"]]
                ax.scatter([row["evaluation"] for row in feasible], [row[key] for row in feasible],
                           color=color, s=6, alpha=.75, linewidths=0)
        ax.set(xlabel="MMA evaluation (cumulative within each run)", ylabel=title)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", color="#dce6ed", linewidth=.8)
        ax.legend(frameon=False, ncols=2)
        with (root / filename).open("xb") as stream:
            fig.savefig(stream, format="png", dpi=150, facecolor="white")
        plt.close(fig)
        images.append(filename)
    return images, limits


def export_comparison(output, cases, *, fitted_spacing_mm=None, write_pngs=True):
    """Replay archived endpoints without optimization and consolidate matched fitted metrics."""
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run this original-FEM export on ASU Sol")
    from shapely.geometry import box

    root = Path(output)
    root.mkdir(parents=True)
    spacing = _latest_common_spacing(cases) if fitted_spacing_mm is None else float(fitted_spacing_mm)
    if not np.isfinite(spacing) or spacing <= 0.:
        raise ValueError("Fitted replay spacing must be finite and positive")
    grid = DesignGrid(.5)
    pipeline = DensityPipeline(grid)
    sharp = box(0., 0., 20., 4.).union(box(9.5, 4., 10.5, 7.))
    sharp_geometry = Reconstruction(sharp, (sharp,), (0, 0, 0), ())
    baseline = {"label": "Original sharp T", "checkpoint": None, "setup": {}, "archive": None}
    results, visuals, histories, slugs = [], [], {}, {"original"}
    for index, case in enumerate((baseline, *cases)):
        label = case["label"]
        slug = "original" if index == 0 else re.sub(r"[^a-z0-9]+", "-", label.lower()).strip("-")
        if index and (not slug or slug in slugs):
            raise ValueError(f"Case label does not uniquely identify a directory: {label}")
        slugs.add(slug)
        folder = root / slug
        folder.mkdir()
        archived = case["checkpoint"]
        row = {"label": label, "archive": case["archive"], "stage": case.get("stage"),
               "flow_band": case["setup"].get("flow_band"), "immersion_operator": case["setup"].get("operator"),
               "archived_kkt_passed": None if archived is None else archived.get("kkt", {}).get("verified"),
               "archived_fitted_screen_passed": None if archived is None else archived.get("fitted_passed"),
               "immersed_flow_pct": None if archived is None else 100 * archived["immersed_side_share"],
               "immersed_power_uW_per_m": None if archived is None else 1e6 * archived["immersed_power_W_per_m"]}
        results.append(row)
        if index:
            history = _history(case)
            histories[label] = history
            with (folder / "mma_history.csv").open("x", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, ["evaluation", "stage", "side_flow_pct", "power_uW_per_m",
                                                 "immersed_constraints_passed"])
                writer.writeheader()
                writer.writerows(history)
        if index and archived is None:
            row["reason"] = case["reason"]
            continue

        if index == 0:
            geometry = sharp_geometry
            projected = grid.base_raw
        else:
            design = np.asarray(archived["design"], dtype=float)
            projected = pipeline.evaluate(design, **case["stage"]).projected
            if archived.get("projected_density") is not None:
                np.testing.assert_allclose(projected, archived["projected_density"], rtol=1e-12, atol=1e-12)
            geometry, replay_area_error = _candidate_geometry(grid, projected, archived)
            row["replay_geometry_symmetric_difference_mm2"] = replay_area_error
        geometry_record = write_reconstruction(geometry, folder / "geometry.json")
        with (folder / "density.npz").open("xb") as stream:
            np.savez_compressed(stream, x_edges_mm=grid.x_edges_mm, y_edges_mm=grid.y_edges_mm,
                                projected_density=projected.reshape(grid.shape))
        if archived is not None and archived.get("geometry", {}).get("geometry_sha256"):
            row["geometry_matches_archive"] = geometry_record["geometry_sha256"] == archived["geometry"]["geometry_sha256"]
        fitted, fields = _fitted_replay(case, geometry, spacing, folder)
        row.update(fitted_replay=fitted, fitted_side_flow_pct=None if "side_share" not in fitted else 100 * fitted["side_share"],
                   fitted_power_uW_per_m=None if "mechanical_power_W_per_m" not in fitted else 1e6 * fitted["mechanical_power_W_per_m"],
                   reconstructed_area_mm2=geometry.area_mm2, connected_ports=geometry.port_component is not None)
        band = row["flow_band"]
        row["fitted_replay_in_flow_band"] = None if band is None or "side_share" not in fitted else (
            band[0] <= fitted["side_share"] <= band[1])
        if index == 0 and (fields is None or not fitted["state_passed"]):
            raise RuntimeError(f"The exact sharp fitted baseline replay did not pass: {fitted}")
        baseline_power = results[0]["fitted_power_uW_per_m"]
        row["fitted_power_change_vs_sharp_pct"] = None if row["fitted_power_uW_per_m"] is None else (
            100. * (row["fitted_power_uW_per_m"] / baseline_power - 1.))
        visuals.append({"folder": folder, "geometry": geometry, "density": projected, "fields": fields})

    columns = ["label", "connected_ports", "reconstructed_area_mm2", "immersed_flow_pct", "immersed_power_uW_per_m",
               "fitted_side_flow_pct", "fitted_power_uW_per_m", "fitted_power_change_vs_sharp_pct",
               "fitted_replay_in_flow_band", "archived_fitted_screen_passed", "archived_kkt_passed", "archive"]
    with (root / "metrics.csv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    report = {"sol_job_id": os.environ["SLURM_JOB_ID"], "baffle": False, "fitted_spacing_mm": spacing,
              "fitted_physics": "Original Taylor-Hood Navier-Stokes FEM, Re=50", "cases": results,
              "image_size_pixels": [1800, 720], "image_extent_mm": VIEW_MM,
              "history_scope": "Archived unique MMA evaluations, stage-separated lines and immersed-feasible dots"}
    if write_pngs:
        try:
            report["images"], report["shared_color_limits"] = _graphics(root, grid, visuals, histories)
        except Exception as error:
            report["image_error"] = f"{type(error).__name__}: {error}"
    _write(root / "comparison.json", report)
    return report


def _mapping(items):
    result = {}
    for item in items:
        label, separator, path = item.partition("=")
        if not separator or not label or not path or label in result:
            raise ValueError("Supply distinct LABEL=PATH pairs")
        result[label] = Path(path)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--lane", action="append", required=True, help="LABEL=LANE_DIRECTORY with setup and stage reports")
    parser.add_argument("--candidate", action="append", default=[], help="Optional LABEL=CHECKPOINT_JSON override")
    parser.add_argument("--fitted-spacing-mm", type=float, help="Default is the finest archived common mesh spacing")
    parser.add_argument("--no-pngs", action="store_true")
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical export in a Sol Slurm job")
    lanes, checkpoints = _mapping(args.lane), _mapping(args.candidate)
    if not checkpoints.keys() <= lanes.keys():
        parser.error("Every candidate override needs a corresponding lane")
    cases = [_load_case(label, lane, checkpoints.get(label)) for label, lane in lanes.items()]
    report = export_comparison(args.output, cases, fitted_spacing_mm=args.fitted_spacing_mm,
                               write_pngs=not args.no_pngs)
    print(json.dumps(report, indent=2, allow_nan=False))
    if report.get("image_error"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
