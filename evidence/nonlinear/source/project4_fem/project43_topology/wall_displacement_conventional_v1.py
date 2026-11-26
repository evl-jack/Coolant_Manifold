"""Boundary-fitted two-rectangle wall-displacement reference for the original T."""

import argparse
from functools import lru_cache
import json
import os
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.optimize import brentq

from .conventional_baseline import MATCH_TOLERANCE_PP, _candidate_checkpoint, _field_archive, _write
from .full_density_continuation import FITTED_POWER_RELATIVE_CHANGE, FITTED_SPACINGS_MM, FITTED_SPLIT_CHANGE_PP
from .reconstruction import Reconstruction, _PORTS_MM, solve_fitted, write_reconstruction


BRANCH_WIDTH_MM = (1.0, 3.0)
UPPER_WALL_Y_MM = (3.0, 5.0)
AREA_MATCH_TOLERANCE_MM2 = 1e-8
BRACKET_SAMPLES = 9


def geometry(branch_width_mm, upper_wall_y_mm):
    """Move only the stated main roof and side-wall segments."""
    from shapely.geometry import LineString, box
    from shapely.ops import unary_union

    b, w = float(branch_width_mm), float(upper_wall_y_mm)
    if not BRANCH_WIDTH_MM[0] <= b <= BRANCH_WIDTH_MM[1]:
        raise ValueError("The branch width must be between 1 and 3 mm")
    if not UPPER_WALL_Y_MM[0] <= w <= UPPER_WALL_Y_MM[1]:
        raise ValueError("The moved main-channel upper wall must be between y=3 and y=5 mm")
    shape = unary_union((box(0.0, 0.0, 12.0, 4.0), box(12.0, 0.0, 15.5, w),
                         box(15.5, 0.0, 20.0, 4.0), box(10.0-b/2.0, 4.0, 10.0+b/2.0, 5.5),
                         box(9.5, 5.5, 10.5, 7.0)))
    if not shape.is_valid or shape.geom_type != "Polygon" or len(shape.interiors):
        raise ValueError("The displaced geometry must be a single valid fluid polygon without holes")
    np.testing.assert_allclose(shape.area, 83.0 + 1.5*(b-1.0) + 3.5*(w-4.0), rtol=0.0, atol=1e-12)
    if not all(shape.boundary.covers(LineString(port)) for port in _PORTS_MM):
        raise ValueError("The displaced geometry must preserve all exact external ports")
    return Reconstruction(shape, (shape,), (0, 0, 0), ())


def _mesh_study(output, b, w, spacings):
    output.mkdir()
    selected = geometry(b, w)
    shape = write_reconstruction(selected, output / "geometry.json")
    rows = []
    for h in spacings:
        start = perf_counter()
        try:
            row, fem, state = solve_fitted(selected, h)
        except Exception as error:
            row = {"spacing_mm": h, "state_passed": False, "error": f"{type(error).__name__}: {error}"}
        else:
            if h == spacings[-1] and row["state_passed"]:
                _field_archive(output / "finest_fitted_fields.npz", fem, state)
        row = {**row, "wall_seconds": perf_counter() - start}
        rows.append(row)
        _write(output / f"fitted_{h:g}_mm.json", row)

    changes = []
    for coarse, fine in zip(rows[:-1], rows[1:]):
        if not coarse["state_passed"] or not fine["state_passed"] or fine["mechanical_power_W_per_m"] <= 0.:
            changes.append({"passed": False, "reason": "A fitted state or its positive power check failed"})
            continue
        split = 100.0 * abs(fine["side_share"] - coarse["side_share"])
        relative_power = abs(coarse["mechanical_power_W_per_m"] / fine["mechanical_power_W_per_m"] - 1.0)
        changes.append({"split_change_pp": split, "relative_power_change": relative_power,
                        "passed": bool(split <= FITTED_SPLIT_CHANGE_PP
                                       and relative_power <= FITTED_POWER_RELATIVE_CHANGE)})
    uncertainty = max(change["split_change_pp"] for change in changes) if all(
        "split_change_pp" in change for change in changes) else None
    report = {"branch_width_mm": b, "upper_wall_y_mm": w, "geometry": shape,
              "physical_fluid_area_mm2": selected.area_mm2, "levels": rows, "successive_changes": changes,
              "split_uncertainty_pp": uncertainty,
              "fitted_mesh_study_passed": all(row["state_passed"] for row in rows)
              and all(change["passed"] for change in changes)}
    _write(output / "mesh_study.json", report)
    return report


def _matched_candidate(output, candidate_path, spacings):
    output.mkdir()
    candidate = _candidate_checkpoint(candidate_path, spacings[-1])
    target_area = float(candidate["physical_fluid_area_mm2"])
    target_split = float(candidate["side_share"])
    if not np.isfinite([target_area, target_split]).all() or not 0.0 <= target_split <= 1.0:
        raise ValueError("The candidate fitted area and split must be finite physical values")

    # Invert each bound on w into a bound on b and intersect with 1 <= b <= 3.
    lower = max(BRANCH_WIDTH_MM[0], 1.0 + (target_area-86.5)/1.5)
    upper = min(BRANCH_WIDTH_MM[1], 1.0 + (target_area-79.5)/1.5)
    report = {"candidate_fitted": candidate, "target_area_mm2": target_area,
              "target_side_share": target_split, "search_area_match_tolerance_mm2": AREA_MATCH_TOLERANCE_MM2,
              "match_tolerance_pp": MATCH_TOLERANCE_PP, "matched": False,
              "qualified_comparison": False}
    if lower > upper:
        return {**report, "reason": "The target area is outside the [79.5, 89.5] mm² family"}
    report["feasible_branch_width_interval_mm"] = [lower, upper]
    trial_dir = output / "search_trials"
    trial_dir.mkdir()
    evaluations, brackets, failures = [], [], []

    def roof(b):
        value = 4.0 + (target_area-83.0-1.5*(b-1.0))/3.5
        # The analytic interval bounds may yield one ulp of roundoff at a wall limit.
        if not UPPER_WALL_Y_MM[0]-1e-12 <= value <= UPPER_WALL_Y_MM[1]+1e-12:
            raise ValueError("The fixed-area branch width produces an infeasible roof position")
        return float(np.clip(value, *UPPER_WALL_Y_MM))

    @lru_cache(maxsize=None)
    def trial(width):
        start = perf_counter()
        b, w = float(width), roof(width)
        shape = geometry(b, w)
        np.testing.assert_allclose(shape.area_mm2, target_area, rtol=0.0, atol=AREA_MATCH_TOLERANCE_MM2)
        try:
            state, _, _ = solve_fitted(shape, spacings[-1])
        except Exception as error:
            state = {"state_passed": False, "error": f"{type(error).__name__}: {error}"}
        residual = state["side_share"] - target_split if state["state_passed"] else None
        row = {"branch_width_mm": b, "upper_wall_y_mm": w, "area_mm2": shape.area_mm2,
               "side_share_residual": residual, "wall_seconds": perf_counter() - start, **state}
        evaluations.append(row)
        _write(trial_dir / f"trial_{len(evaluations):03d}.json", row)
        return residual

    def bracket_function(width):
        residual = trial(float(width))
        if residual is None:
            raise RuntimeError(f"A fitted solve failed inside the bracket at b={width:.9g} mm")
        return residual

    last = None
    selected = None
    for b in np.unique(np.linspace(lower, upper, BRACKET_SAMPLES)):
        b = float(b)
        residual = trial(b)
        if residual is None:
            last = None
        else:
            if abs(residual) <= 1e-10 and selected is None:
                selected = b
            if last is not None and residual * last[1] < 0.0:
                brackets.append((last[0], b))
            last = (b, residual)

    used_bracket = None
    for bracket in brackets:
        if selected is not None:
            break
        try:
            selected = float(brentq(bracket_function, *bracket, xtol=1e-5, rtol=1e-9, maxiter=40))
        except (ValueError, RuntimeError) as error:
            failures.append({"bracket_mm": bracket, "error": f"{type(error).__name__}: {error}"})
        else:
            used_bracket = bracket

    report.update(search_history=evaluations, candidate_brackets_mm=brackets, bracket_failures=failures,
                  used_bracket_mm=used_bracket)
    if selected is None:
        passing = [row for row in evaluations if row["state_passed"]]
        if not passing:
            return {**report, "reason": "No fitted state passed the bounded sweep"}
        nearest = min(passing, key=lambda row: abs(row["side_share_residual"]))
        report["closest_valid_sample"] = nearest
        if 100.0 * abs(nearest["side_share_residual"]) > MATCH_TOLERANCE_PP:
            return {**report, "reason": "No split match was found in the sampled fixed-area interval"}
        selected = nearest["branch_width_mm"]

    mesh = _mesh_study(output / "matched_geometry", selected, roof(selected), spacings)
    report["conventional"] = mesh
    final = mesh["levels"][-1]
    if not final["state_passed"]:
        return {**report, "reason": "The final independent fitted replay failed"}
    split_error = 100.0 * (final["side_share"] - target_split)
    area_error = mesh["physical_fluid_area_mm2"] - target_area
    matched = abs(split_error) <= MATCH_TOLERANCE_PP and abs(area_error) <= AREA_MATCH_TOLERANCE_MM2
    qualified = matched and mesh["fitted_mesh_study_passed"] and candidate["physical_screen_passed"]
    report.update(matched=matched, split_error_pp=split_error, area_error_mm2=area_error,
                  qualified_comparison=qualified)
    if qualified:
        report["candidate_power_change_vs_conventional_pct"] = 100.0 * (
            candidate["mechanical_power_W_per_m"] / final["mechanical_power_W_per_m"] - 1.0)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="New archive directory")
    parser.add_argument("--candidate", type=Path, help="Candidate report with independently fitted levels")
    parser.add_argument("--fit-spacings-mm", nargs=3, type=float, default=FITTED_SPACINGS_MM)
    args = parser.parse_args()
    job_id = os.environ.get("SLURM_JOB_ID")
    if not job_id:
        parser.error("Run this fitted comparison in a Sol Slurm job")
    spacings = tuple(args.fit_spacings_mm)
    if (any(not np.isfinite(h) or h <= 0.0 for h in spacings)
            or any(h <= finer for h, finer in zip(spacings[:-1], spacings[1:]))):
        parser.error("Fitted spacings must be positive and strictly decreasing")

    args.output.mkdir(parents=True)
    report = {"sol_job_id": job_id, "fit_spacings_mm": spacings,
              "split_mesh_change_limit_pp": FITTED_SPLIT_CHANGE_PP,
              "power_mesh_change_limit_relative": FITTED_POWER_RELATIVE_CHANGE}
    try:
        if args.candidate is None:
            initial = _mesh_study(args.output / "b3_w3", 3.0, 3.0, spacings)
            report["b3_w3"] = initial
            passed = initial["fitted_mesh_study_passed"]
        else:
            matched = _matched_candidate(args.output / "area_and_split_match", args.candidate, spacings)
            report["candidate_match"] = matched
            passed = matched["qualified_comparison"]
    except Exception as error:
        _write(args.output / "failure.json", {**report, "error": f"{type(error).__name__}: {error}"})
        raise
    _write(args.output / "wall_displacement_conventional_v1.json", {**report, "passed": passed})
    raise SystemExit(0 if passed else 2)


if __name__ == "__main__":
    main()
