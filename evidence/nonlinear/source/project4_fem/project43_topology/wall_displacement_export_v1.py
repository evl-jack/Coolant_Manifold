"""Export archived wall-displacement PNGs and matched fitted-wall evidence."""

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path

import numpy as np

from .baffle_free_export import (
    DEFAULT_FITTED_SPACING_MM, VIEW_MM, _graphics, _read, _state_archive, _write,
)
from .density import DesignGrid
from .reconstruction import Reconstruction, solve_fitted, write_reconstruction
from .wall_displacement_density_v1 import wall_displacement_mask
from .wall_displacement_geometry_v1 import reconstruct_clipped, sharp_geometry


HISTORY_COLUMNS = (
    "evaluation", "stage", "stage_evaluation", "objective", "power_W_per_m",
    "main_flow_pct", "side_flow_pct", "lower_inequality", "upper_inequality",
    "optimizer_feasible",
)


def _checkpoint(lane, summary):
    final = summary.get("final_candidate")
    if final is not None:
        path = lane / final / "checkpoint.json"
        archived = _read(path)
        if not summary["completed_schedule"] or not archived.get("qualified_optimization_terminal"):
            raise ValueError("The reported final checkpoint did not complete the optimization schedule")
        return path, archived, "qualified_optimization_endpoint"
    for stage in sorted(lane.glob("stage_[0-9][0-9]"), reverse=True):
        for label in ("terminal", "best_feasible_trial", "restored_start"):
            path = stage / label / "checkpoint.json"
            if path.is_file():
                archived = _read(path)
                if "projected_density" in archived:
                    return path, archived, f"diagnostic_{label}"
    return None, None, "no_archived_candidate"


def _live_stages(lane):
    stages = {}
    for path in sorted(lane.glob("stage_[0-9][0-9]/summary.json")):
        payload = path.read_text(encoding="utf-8")
        # The runner writes its JSON reports directly. A concurrent read may
        # see the newest summary before its final newline has been written.
        if not payload.endswith("\n"):
            continue
        report = json.loads(payload)
        if report["completed"]:
            stages[path.parent] = report["results"]
    return stages


def _live_checkpoint(stages):
    for stage, summary in reversed(stages.items()):
        path = stage / "terminal" / "checkpoint.json"
        if not path.is_file():
            continue
        archived = _read(path)
        if "projected_density" not in archived or not (path.parent / archived["field_archive"]).is_file():
            continue
        if archived["stage"] != summary["stage"]:
            raise ValueError(f"The terminal checkpoint disagrees with its stage summary: {path}")
        qualified = archived["qualified_optimization_terminal"]
        selection = "snapshot_stage_terminal" if qualified else "snapshot_diagnostic_terminal"
        return path, archived, selection
    return None, None, "snapshot_no_archived_terminal"


def _history(lane, setup, *, stages=None):
    low, high = setup["flow_band"]
    scale = setup["flow_constraint_scale"]
    reference = setup["normalization_W_per_m"]
    rows = []
    for path in sorted(lane.glob("stage_[0-9][0-9]/mma.json")):
        if stages is not None and path.parent not in stages:
            continue
        outcome = _read(path)
        stage = int(path.parent.name.split("_")[1])
        for index, row in enumerate(outcome["evaluation_history"], start=1):
            lower, upper = row["inequalities"]
            fraction = low - scale * lower
            if not np.isclose(fraction, high + scale * upper, rtol=0., atol=1e-10):
                raise ValueError(f"Inconsistent flow inequalities in {path}")
            rows.append({
                "evaluation": len(rows) + 1, "stage": stage, "stage_evaluation": index,
                "objective": row["objective"], "power_W_per_m": reference * row["objective"],
                "main_flow_pct": 100. * (1. - fraction), "side_flow_pct": 100. * fraction,
                "lower_inequality": lower, "upper_inequality": upper,
                "optimizer_feasible": row["feasible"],
            })
    return rows


def _checked_candidate(grid, free, archive):
    raw = np.asarray(archive["raw_density"], dtype=float)
    projected = np.asarray(archive["projected_density"], dtype=float)
    design = np.asarray(archive["design"], dtype=float)
    for field in (raw, projected):
        if field.shape != (grid.size,) or not np.isfinite(field).all():
            raise ValueError("An archived density does not cover the full 960-cell grid")
        if np.any((field < -1e-12) | (field > 1. + 1e-12)):
            raise ValueError("An archived density is outside [0, 1]")
    if not np.array_equal(raw[free], design):
        raise ValueError("The archived compact controls disagree with the raw field")
    if not np.array_equal(raw[~free], grid.base_raw[~free]):
        raise ValueError("The archived raw field changed an anchored cell")
    np.testing.assert_allclose(projected[~free], grid.base_raw[~free], atol=1e-12, rtol=0.)

    reference = archive.get("geometry")
    if reference is None:
        return raw, projected, None, None
    from shapely.geometry import shape

    geojson = reference["geometry"]
    if sha256(json.dumps(geojson, sort_keys=True).encode()).hexdigest() != reference["geometry_sha256"]:
        raise ValueError("The archived clipped geometry digest did not verify")
    selected = shape(geojson)
    replay = reconstruct_clipped(grid, projected, free)
    difference = float(selected.symmetric_difference(replay.geometry).area)
    if difference > 1e-10 or tuple(reference["port_components"]) != replay.port_components:
        raise ValueError("The archived geometry disagrees with the projected-density clipping replay")
    polys = [selected] if selected.geom_type == "Polygon" else selected.geoms
    ordered = tuple(sorted(polys, key=lambda poly: (*poly.bounds, poly.area)))
    geometry = Reconstruction(selected, ordered, replay.port_components, replay.contour_segments_mm)
    return raw, projected, geometry, difference


def _fitted(geometry, spacing, folder, *, disabled):
    if disabled:
        return {"attempted": False, "reason": "Fitted replay was disabled"}, None
    if geometry is None or geometry.port_component is None:
        return {"attempted": False, "reason": "No archived geometry connects all three ports"}, None
    try:
        metrics, fem, state = solve_fitted(geometry, spacing)
    except (RuntimeError, ValueError, AssertionError, np.linalg.LinAlgError) as error:
        return {"attempted": True, "error": f"{type(error).__name__}: {error}"}, None
    if not metrics["state_passed"]:
        return dict(metrics, attempted=True), None
    return dict(metrics, attempted=True), _state_archive(folder / "fitted_fields.npz", fem, state)


def _cell_image(path, grid, field, *, mask=False):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from shapely.geometry import box

    fig = plt.figure(figsize=(12, 4.8), dpi=150, facecolor="white")
    try:
        ax = fig.add_axes((0., 0., 1., 1.))
        ax.set(aspect="equal", xlim=VIEW_MM[:2], ylim=VIEW_MM[2:])
        ax.set_axis_off()
        cmap = plt.get_cmap("viridis")
        if mask:
            boundary = sharp_geometry().exterior
            x, y = boundary.xy
            ax.plot(x, y, color="#17415c", linewidth=1., zorder=3)
            cmap = ListedColormap(["#edf2f4", "#da7735"])
        else:
            for collar in (box(0, 0, 6, 4), box(16, 0, 20, 4), box(9.5, 6, 10.5, 7)):
                x, y = collar.exterior.xy
                ax.fill(x, y, color=cmap(1.), linewidth=0, antialiased=False)
        ax.pcolormesh(grid.x_edges_mm, grid.y_edges_mm, field.reshape(grid.shape),
                      cmap=cmap, vmin=0, vmax=1, shading="flat", antialiased=False)
        with path.open("xb") as stream:
            fig.savefig(stream, format="png", dpi=150, facecolor="white")
    finally:
        plt.close(fig)


def _history_images(output, rows, band):
    import matplotlib.pyplot as plt

    images = []
    for name, series, ylabel in (
        ("objective_history.png", (("objective", "Normalized power"),), "Pumping power / sharp reference"),
        ("flow_constraints_history.png", (("side_flow_pct", "Side outlet"), ("main_flow_pct", "Main outlet")),
         "Share of inlet flow (%)"),
    ):
        fig, ax = plt.subplots(figsize=(12, 4.8), dpi=150, layout="constrained")
        try:
            if name.startswith("flow"):
                ax.axhspan(100. * band[0], 100. * band[1], color="#e2f1e8", label="Allowed flow band")
            else:
                ax.axhline(1., color="#77624d", linewidth=1., linestyle="--", label="Sharp fitted reference")
            for stage_index, stage in enumerate(sorted({row["stage"] for row in rows})):
                stage_rows = [row for row in rows if row["stage"] == stage]
                x = [row["evaluation"] for row in stage_rows]
                if stage_index:
                    ax.axvline(x[0] - .5, color="#d0d9df", linewidth=.7)
                for (key, label), color in zip(series, ("#1b577a", "#b34f25")):
                    ax.plot(x, [row[key] for row in stage_rows], color=color, linewidth=1.2,
                            marker=".", markersize=2.5, label=label if stage_index == 0 else None)
            ax.set(xlabel="Unique MMA evaluation", ylabel=ylabel)
            ax.spines[["top", "right"]].set_visible(False)
            ax.grid(axis="y", color="#dce6ed", linewidth=.8)
            ax.legend(frameon=False, ncols=3)
            with (output / name).open("xb") as stream:
                fig.savefig(stream, format="png", dpi=150, facecolor="white")
        finally:
            plt.close(fig)
        images.append(name)
    return images


def export(output, lane=None, *, spacing=DEFAULT_FITTED_SPACING_MM, skip_fitted=False,
           physical_audit=None, live_snapshot=False):
    """Export exact archived densities and geometry with optional matched-spacing FEM."""
    job = os.environ.get("SLURM_JOB_ID")
    if not job:
        raise RuntimeError("Run this export as a Sol Slurm job")
    if not np.isfinite(spacing) or spacing <= 0.:
        raise ValueError("Fitted spacing must be finite and positive")
    if live_snapshot and (lane is None or physical_audit is not None):
        raise ValueError("Live snapshots require an optimization lane and no final physical audit")
    grid = DesignGrid(.25)
    free = wall_displacement_mask(grid)
    summary, archived, checkpoint, history, stages = {}, None, None, [], {}
    selection = "no_optimization_archive"
    band = (.10, .90)
    if lane is not None:
        setup = _read(lane / "setup.json")
        if (setup["design_cells"] != grid.size or setup["projection_threshold"] != .5
                or setup["design_spacing_mm"] != grid.spacing_mm
                or not np.array_equal(setup["free_cells"], free)
                or not np.array_equal(setup["original_raw_sharp_t_design"], grid.base_raw)):
            raise ValueError("The archived setup does not match the v1 design grid and initialization")
        band = tuple(setup["flow_band"])
        if band != (.10, .90):
            raise ValueError("The archived flow band differs from the v1 limits")
        if live_snapshot:
            stages = _live_stages(lane)
            summary = next(reversed(stages.values()), {})
            checkpoint, archived, selection = _live_checkpoint(stages)
            history = _history(lane, setup, stages=stages)
        else:
            summary_path = lane / "summary.json"
            if summary_path.is_file():
                wrapped = json.loads(summary_path.read_text(encoding="utf-8"))
                if wrapped["completed"]:
                    summary = wrapped["results"]
            checkpoint, archived, selection = _checkpoint(lane, summary)
            history = _history(lane, setup)

    output.mkdir(parents=True)
    original = output / "original"
    original.mkdir()
    sharp = sharp_geometry()
    sharp = Reconstruction(sharp, (sharp,), (0, 0, 0), ())
    write_reconstruction(sharp, original / "geometry.json")
    fitted, fields = _fitted(sharp, spacing, original, disabled=skip_fitted)
    visuals = [dict(folder=original, geometry=sharp, density=grid.base_raw, fields=fields)]
    images = ["editable_mask.png"]
    _cell_image(output / images[0], grid, free, mask=True)
    materials = {
        "x_edges_mm": grid.x_edges_mm, "y_edges_mm": grid.y_edges_mm,
        "editable_mask": free.reshape(grid.shape), "original_binary_density": grid.base_raw.reshape(grid.shape),
    }
    report = {
        "sol_job_id": job, "selection": selection, "lane": str(lane) if lane else None,
        "checkpoint": str(checkpoint) if checkpoint else None, "optimizer_stop_reason": summary.get("stop_reason"),
        "sharp_fitted": fitted, "matched_nominal_fitted_spacing_mm": spacing, "flow_band": band,
        "image_size_pixels": [1800, 720], "image_extent_mm": VIEW_MM, "design_cells": grid.size,
        "editable_cells": int(free.sum()), "projection_threshold": .5,
        "history_scope": "Archived unique MMA evaluations. Negative inequalities satisfy the constraints.",
        "fitted_scope": "Single-spacing original Taylor-Hood Navier-Stokes at Re=50 for visualization",
    }
    if live_snapshot:
        report["snapshot_completed_stages"] = [stage.name for stage in stages]
    if physical_audit is not None:
        physical = _read(physical_audit)
        if physical.get("attempted"):
            if (selection != "qualified_optimization_endpoint" or
                    Path(physical["checkpoint"]).resolve() != checkpoint.resolve()):
                raise ValueError("The independent physical audit refers to a different checkpoint")
            if physical["geometry"]["geometry_sha256"] != (archived.get("geometry") or {}).get("geometry_sha256"):
                raise ValueError("The independently fitted polygon differs from the optimizer archive")
        report["independent_physical_audit"] = {
            "archive": str(physical_audit), "attempted": physical["attempted"], "reason": physical.get("reason"),
            "physical_screen_passed": physical.get("physical_screen_passed"),
            "smoothness_passed": physical.get("smoothness_passed"),
        }
    if archived is not None:
        candidate = output / selection
        candidate.mkdir()
        raw, projected, geometry, difference = _checked_candidate(grid, free, archived)
        materials.update(archived_raw_density=raw.reshape(grid.shape),
                         archived_projected_density=projected.reshape(grid.shape))
        raw_image = candidate / "raw_density.png"
        _cell_image(raw_image, grid, raw)
        images.append(raw_image.relative_to(output).as_posix())
        candidate_fitted, fields = _fitted(geometry, spacing, candidate, disabled=skip_fitted)
        if geometry is None:
            projected_image = candidate / "projected_density.png"
            _cell_image(projected_image, grid, projected)
            images.append(projected_image.relative_to(output).as_posix())
        else:
            written = write_reconstruction(geometry, candidate / "geometry.json")
            if written["geometry_sha256"] != archived["geometry"]["geometry_sha256"]:
                raise ValueError("The exported clipped geometry digest differs from the archive")
            visuals.append(dict(folder=candidate, geometry=geometry, density=projected, fields=fields))
        report["candidate"] = {
            "stage": archived["stage"], "optimization_terminal_qualified": archived.get("qualified_optimization_terminal"),
            "wall_geometry_audit": archived.get("wall_geometry_audit"), "geometry_area_difference_mm2": difference,
            "immersed_side_share": archived["side_fraction"], "immersed_power_W_per_m": archived["power_w_per_m"],
            "fitted_replay": candidate_fitted,
            "fitted_flow_band_passed": None if not candidate_fitted.get("state_passed") else
            bool(band[0] <= candidate_fitted["side_share"] <= band[1]),
        }
        if fitted.get("state_passed") and candidate_fitted.get("state_passed"):
            power_ratio = candidate_fitted["mechanical_power_W_per_m"] / fitted["mechanical_power_W_per_m"]
            report["candidate"]["fitted_power_change_vs_sharp_pct"] = 100. * (power_ratio - 1.)
    with (output / "density_fields.npz").open("xb") as stream:
        np.savez_compressed(stream, **materials)
    with (output / "mma_history.csv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, HISTORY_COLUMNS)
        writer.writeheader()
        writer.writerows(history)
    panels, report["shared_color_limits"] = _graphics(output, grid, visuals, {})
    # The shared renderer writes density.png for all panels; expose the actual field names.
    for index, name in enumerate(panels):
        if Path(name).name == "density.png":
            path = output / name
            renamed = path.with_name("binary_density.png" if path.parent == original else "projected_density.png")
            path.rename(renamed)
            panels[index] = renamed.relative_to(output).as_posix()
    images.extend(Path(path).as_posix() for path in panels)
    if history:
        images.extend(_history_images(output, history, band))
    report["images"] = images
    _write(output / "export.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="New export directory")
    parser.add_argument("--lane", type=Path, help="Wall-displacement optimization archive")
    parser.add_argument("--fitted-spacing-mm", type=float, default=DEFAULT_FITTED_SPACING_MM)
    parser.add_argument("--skip-fitted", action="store_true", help="Export archived geometry and histories only")
    parser.add_argument("--physical-audit", type=Path, help="Optional independent physical_audit.json")
    parser.add_argument("--live-snapshot", action="store_true",
                        help="Read terminal states and MMA history from completed stage summaries")
    args = parser.parse_args()
    if args.physical_audit is not None and args.lane is None:
        parser.error("--physical-audit requires an optimization archive in --lane")
    if args.live_snapshot and (args.lane is None or args.physical_audit is not None):
        parser.error("--live-snapshot requires --lane and cannot be used with --physical-audit")
    export(args.output, args.lane, spacing=args.fitted_spacing_mm,
           skip_fitted=args.skip_fitted, physical_audit=args.physical_audit, live_snapshot=args.live_snapshot)


if __name__ == "__main__":
    main()
