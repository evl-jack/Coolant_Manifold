"""Fit an independently meshed, centered main-channel contraction to a Model B split."""

import argparse
from functools import lru_cache
import json
import os
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import brentq

from .full_density_continuation import FITTED_POWER_RELATIVE_CHANGE, FITTED_SPLIT_CHANGE_PP
from .reconstruction import Reconstruction, _PORTS_MM, solve_fitted, write_reconstruction


SIDE_BAND = (0.10, 0.90)
AREA_CAP_MM2 = 83.0
MATCH_TOLERANCE_PP = 0.02
BRACKET_WIDTHS_MM = (4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.60, 0.30, 0.15)


def geometry(throat_width_mm):
    """Preserve the original sharp side leg and all exact port boundaries."""
    from shapely.geometry import LineString, Polygon, box

    if not np.isfinite(throat_width_mm) or not 0.15 <= throat_width_mm <= 4.0:
        raise ValueError("Throat width must lie in [0.15, 4] mm")
    h = float(throat_width_mm)
    x = (0.0, 11.5, 13.0, 14.0, 16.0, 20.0)
    bottom = (0.0, 0.0, 2.0-h/2, 2.0-h/2, 0.0, 0.0)
    lower = list(zip(x, bottom))
    upper = [(xi, 4.0-yi) for xi, yi in reversed(lower)]
    shape = Polygon(lower + upper).union(box(9.5, 4.0, 10.5, 7.0))
    if not shape.is_valid or shape.geom_type != "Polygon":
        raise ValueError("The conventional design must have a single valid fluid polygon")
    np.testing.assert_allclose(shape.area, 83.0-2.75*(4.0-h), rtol=1e-12, atol=1e-12)
    if not all(shape.boundary.covers(LineString(port)) for port in _PORTS_MM):
        raise ValueError("The conventional geometry changed an external port")
    return Reconstruction(shape, (shape,), (0, 0, 0), ())


def _write(path, value):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def _field_archive(path, fem, state):
    from .domain import REFERENCE_HEIGHT_MM
    from .evaluator import HEIGHT_M, VELOCITY_SCALE, VISCOSITY

    with path.open("xb") as stream:
        np.savez_compressed(stream, vertices_mm=fem.points * REFERENCE_HEIGHT_MM, triangles=fem.triangles,
                            velocity_vertices_mm=fem.velocity_points * REFERENCE_HEIGHT_MM,
                            velocity_dofs=fem.velocity_dofs,
                            velocity_m_s=VELOCITY_SCALE * state.values[:2*fem.n_u].reshape((2, fem.n_u)).T,
                            pressure_pa=VISCOSITY*VELOCITY_SCALE/HEIGHT_M * state.values[2*fem.n_u:])


def _geometry_figures(output, selected):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    figures = []
    for name, polygon in (("original_sharp_t", geometry(4.).port_component),
                          ("conventional_geometry", selected.port_component)):
        fig = plt.figure(figsize=(12., 4.8), dpi=150, facecolor="white")
        try:
            ax = fig.add_axes((0., 0., 1., 1.))
            ax.set(aspect="equal", xlim=(-.5, 20.5), ylim=(-.7, 7.7))
            ax.set_axis_off()
            ax.add_patch(Polygon(polygon.exterior.coords, closed=True, facecolor="#17415c", edgecolor="none"))
            file = f"{name}.png"
            with (output / file).open("xb") as stream:
                fig.savefig(stream, format="png", dpi=150, facecolor="white")
            figures.append(file)
        finally:
            plt.close(fig)
    return figures


def _candidate_checkpoint(path, spacing):
    archived = json.loads(path.read_text(encoding="utf-8"))
    candidate = archived.get("results", archived)
    fitted = candidate.get("fitted", candidate.get("fitted_screen", {}))
    rows = fitted.get("levels", [])
    matching = [row for row in rows if np.isclose(row["spacing_mm"], spacing, atol=1e-12, rtol=0.)]
    if not matching or not matching[-1]["state_passed"]:
        raise ValueError("The target candidate needs an archived passing fitted state at the matched spacing")
    row = dict(matching[-1])
    row["archive"] = str(path)
    row["physical_screen_passed"] = fitted.get("fitted_diagnostic_passed", False)
    return row


def run(output, target, spacings, candidate=None):
    """Tune a geometric reference to the same flow split under unchanged traction conditions."""
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run this numerical comparison in a Sol Slurm job")
    spacings = tuple(float(h) for h in spacings)
    if (len(spacings) != 3 or any(not np.isfinite(h) or h <= 0. for h in spacings)
            or any(a <= b for a, b in zip(spacings[:-1], spacings[1:]))):
        raise ValueError("Supply three finite, positive and strictly decreasing fitted spacings")
    if not np.isfinite(target) or not SIDE_BAND[0] <= target <= SIDE_BAND[1]:
        raise ValueError("Target side fraction must lie in Model B's 10–90% band")
    trial_dir = output / "search_trials"
    trial_dir.mkdir()
    evaluations = []

    @lru_cache(maxsize=None)
    def trial(width):
        start = perf_counter()
        width = float(width)
        try:
            record, fem, state = solve_fitted(geometry(width), spacings[-1])
        except Exception as error:
            record = {"error": f"{type(error).__name__}: {error}"}
        entry = {"throat_width_mm": width, "wall_seconds": perf_counter()-start, **record}
        evaluations.append(entry)
        _write(trial_dir / f"trial_{len(evaluations):03d}.json", entry)
        if not entry.get("state_passed", False):
            raise RuntimeError(f"Fitted search failed at {width:.8g} mm: {entry}")
        if width == 4.:
            _field_archive(output / "original_sharp_t_fitted_fields.npz", fem, state)
        return record["side_share"] - target

    bracket = None
    width = None
    previous_width, previous_residual = None, None
    for trial_width in BRACKET_WIDTHS_MM:
        residual = trial(trial_width)
        if abs(residual) <= 1e-10:
            width = trial_width
            break
        if previous_residual is not None and residual * previous_residual < 0:
            bracket = (trial_width, previous_width)
            break
        previous_width, previous_residual = trial_width, residual
    if width is None and bracket is None:
        return {"matched": False, "reason": "The geometric family did not bracket the requested split",
                "target_side_fraction": target, "search_history": evaluations,
                "sol_job_id": os.environ["SLURM_JOB_ID"]}
    if width is None:
        width = float(brentq(trial, *bracket, xtol=1e-5, rtol=1e-9, maxiter=40))

    selected = geometry(width)
    shape = write_reconstruction(selected, output / "conventional_geometry.json")
    reports = []
    for spacing in spacings:
        start = perf_counter()
        row, fem, state = solve_fitted(selected, spacing)
        reports.append(dict(row, wall_seconds=perf_counter()-start))
        _write(output / f"fitted_{spacing:g}_mm.json", reports[-1])
    _field_archive(output / "conventional_fitted_fields.npz", fem, state)
    changes = []
    for a, b in zip(reports[:-1], reports[1:]):
        split_pp = 100 * abs(b["side_share"] - a["side_share"])
        power = b["mechanical_power_W_per_m"]
        relative_power = abs(power - a["mechanical_power_W_per_m"]) / power if power > 0. else None
        changes.append(dict(split_change_pp=split_pp, relative_power_change=relative_power,
                            passed=split_pp <= FITTED_SPLIT_CHANGE_PP and relative_power is not None
                            and relative_power <= FITTED_POWER_RELATIVE_CHANGE))
    final = reports[-1]
    split_error_pp = 100 * (final["side_share"] - target)
    uncertainty = max(change["split_change_pp"] for change in changes) / 100.
    split_interval = (final["side_share"]-uncertainty, final["side_share"]+uncertainty)
    in_band = SIDE_BAND[0] <= split_interval[0] and split_interval[1] <= SIDE_BAND[1]
    match = abs(split_error_pp) <= MATCH_TOLERANCE_PP
    result = {"sol_job_id": os.environ["SLURM_JOB_ID"], "matched": match,
              "target_side_fraction": target, "match_error_pp": split_error_pp,
              "match_tolerance_pp": MATCH_TOLERANCE_PP, "throat_width_mm": width,
              "bracket_mm": bracket, "geometry": shape, "search_history": evaluations,
              "levels": reports, "successive_changes": changes, "side_fraction_interval": split_interval,
              "fitted_fields": {"sharp": "original_sharp_t_fitted_fields.npz",
                                "conventional": "conventional_fitted_fields.npz"},
              "same_mesh_sharp_fitted": evaluations[0],
              "conventional_power_change_vs_sharp_pct": 100. * (
                  final["mechanical_power_W_per_m"] / evaluations[0]["mechanical_power_W_per_m"] - 1.),
              "side_fraction_interval_in_model_b_band": in_band,
              "physical_model_b_passed": bool(selected.area_mm2 <= AREA_CAP_MM2 and in_band
                                             and all(r["state_passed"] for r in reports)
                                             and all(c["passed"] for c in changes))}
    try:
        result["geometry_images"] = _geometry_figures(output, selected)
    except Exception as error:
        result["geometry_image_error"] = f"{type(error).__name__}: {error}"
    if candidate is not None:
        result["candidate_fitted"] = candidate
        result["candidate_power_change_vs_conventional_pct"] = 100 * (
            candidate["mechanical_power_W_per_m"] / final["mechanical_power_W_per_m"] - 1.)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="New archive directory")
    targets = parser.add_mutually_exclusive_group()
    targets.add_argument("--target-side-fraction", type=float, default=0.11)
    targets.add_argument("--candidate", type=Path, help="Model B endpoint report with independently fitted levels")
    parser.add_argument("--fit-spacings-mm", nargs=3, type=float, default=(0.5, 0.25, 0.125))
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical comparison in a Sol Slurm job")
    candidate = _candidate_checkpoint(args.candidate, args.fit_spacings_mm[-1]) if args.candidate else None
    target = candidate["side_share"] if candidate is not None else args.target_side_fraction
    args.output.mkdir(parents=True)
    try:
        report = run(args.output, target, args.fit_spacings_mm, candidate)
    except Exception as error:
        _write(args.output / "failure.json", {"error": f"{type(error).__name__}: {error}"})
        raise
    _write(args.output / "conventional_baseline.json", report)
    print(json.dumps({k: v for k, v in report.items() if k not in ("geometry", "search_history")},
                     indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
