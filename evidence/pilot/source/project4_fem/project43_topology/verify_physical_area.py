"""Independently verify the reconstructed-area upper bound on Sol."""

import argparse
import json
import os
from pathlib import Path
import traceback

import numpy as np

from .density import DesignGrid, DensityPipeline, FIXED_FLUID_AREA_MM2
from .expanded_envelope import ExpandedDesignGrid
from .gate2 import write_report
from .physical_area import ReconstructedAreaUpperBound
from .reconstruction import reconstruct
from .verify_density import DERIVATIVE_TOLERANCE, MINIMUM_TAYLOR_ORDER
from .verify_flow import EQUIVALENCE_TOLERANCE


ARCHIVE_RELATIVE_PATH = (
    "p43-relaxed-flow.rsOvla/original-fem-model-b-2405258/stage_00/split_margin/checkpoint.json"
)
AREA_CAP_MM2 = 83.0
TRANSITION = 0.02
REGULARIZATION = 1e-6


def _geometry_check(grid, area_model, projected):
    approximate, gradient = area_model.evaluate(projected)
    exact = reconstruct(grid, projected).area_mm2
    gap = approximate - exact
    return {
        "bound_mm2": approximate,
        "reconstructed_mm2": exact,
        "bound_minus_reconstructed_mm2": gap,
        "bound_margin_to_cap_mm2": AREA_CAP_MM2 - approximate,
        "reconstructed_margin_to_cap_mm2": AREA_CAP_MM2 - exact,
        "gradient_inf_norm": float(np.max(np.abs(gradient))),
        "passed": bool(np.isfinite(approximate) and np.isfinite(gradient).all()
                       and gap >= -EQUIVALENCE_TOLERANCE),
    }


def _direction_check(response, values, direction, *, at_bound=False):
    base, gradient = response(values)
    analytic = gradient @ direction
    used = direction != 0.0
    if not used.any():
        raise ValueError("A verification direction must be nonzero")
    if at_bound:
        feasible_distance = np.where(direction[used] > 0.0, 1.0 - values[used], values[used])
        feasible_step = float(np.min(feasible_distance / np.abs(direction[used]))) / 2.0
    else:
        feasible_distance = np.minimum(values[used], 1.0 - values[used])
        feasible_step = float(np.min(feasible_distance / np.abs(direction[used])))
    initial_step = min(0.01, 0.8 * feasible_step)
    if initial_step <= 0.0:
        raise ValueError("The perturbation must remain in the density box")
    steps = initial_step / 2.0**np.arange(20)
    samples = []
    for step in steps:
        plus = response(values + step * direction)[0]
        if at_bound:
            ahead = response(values + 2.0 * step * direction)[0]
            derivative = (4.0 * (plus - base) - (ahead - base)) / (2.0 * step)
        else:
            minus = response(values - step * direction)[0]
            derivative = (plus - minus) / (2.0 * step)
        samples.append({
            "step": float(step),
            "fd_area_derivative_mm2": float(derivative),
            "scaled_derivative_error": float(abs(derivative - analytic) / max(abs(analytic), 1.0)),
            "taylor_remainder_mm2": float(abs(plus - base - step * analytic)),
        })

    remainders = [row["taylor_remainder_mm2"] for row in samples]
    floor = 1e-12
    orders = [float(np.log2(broad / narrow)) if min(broad, narrow) > floor else None
              for broad, narrow in zip(remainders[:-1], remainders[1:])]
    affine_to_roundoff = all(remainder <= floor for remainder in remainders)
    order_passed = affine_to_roundoff or any(
        first is not None and second is not None and min(first, second) >= MINIMUM_TAYLOR_ORDER
        for first, second in zip(orders[:-1], orders[1:])
    )
    # Two successive samples near the local asymptotic window are required.
    # Smaller steps remain archived to show when cancellation takes over.
    adequate = [row["scaled_derivative_error"] <= DERIVATIVE_TOLERANCE for row in samples[6:]]
    finite_difference_passed = any(left and right for left, right in zip(adequate[:-1], adequate[1:]))
    return {
        "difference": "second-order inward" if at_bound else "centered",
        "analytic_area_derivative_mm2": float(analytic),
        "samples": samples,
        "taylor_orders": orders,
        "affine_to_roundoff": affine_to_roundoff,
        "finite_difference_passed": finite_difference_passed,
        "taylor_passed": order_passed,
        "passed": order_passed and finite_difference_passed,
    }


def _knot_fixtures(grid):
    saddle = grid.base_raw.reshape(grid.shape).copy()
    checker_shape = np.array([[0.75, 0.25], [0.25, 0.75]])
    saddle[1:3, 1:3] = checker_shape
    saddle_perturbation = np.zeros(grid.shape)
    saddle_perturbation[1, 1] = 1.0
    switch = saddle.copy()
    switch[1:3, 1:3] = np.array([[0.2, 0.8], [0.7, 0.5-TRANSITION]])
    switch_perturbation = np.zeros(grid.shape)
    switch_perturbation[2, 2] = 1.0
    tie = saddle.copy()
    tie[1:3, 1:3] = np.array([[0.5, 0.25], [0.75, 0.85]])
    tie_perturbation = np.zeros(grid.shape)
    tie_perturbation[1, 1] = 1.0
    return {"checkerboard_saddle": (saddle.ravel(), saddle_perturbation.ravel()),
            "transition_start": (switch.ravel(), switch_perturbation.ravel()),
            "exact_threshold": (tie.ravel(), tie_perturbation.ravel())}


def _envelope_check(grid, rng):
    area = ReconstructedAreaUpperBound(grid, transition=TRANSITION, regularization=REGULARIZATION)
    projected_random = rng.uniform(0.12, 0.9, grid.size)
    near_threshold = rng.uniform(0.5 - 1.25 * TRANSITION, 0.5 + 1.25 * TRANSITION, grid.size)
    projected_fixtures = {
        "all_solid": np.zeros(grid.size),
        "all_fluid": np.ones(grid.size),
        "all_threshold_ties": np.full(grid.size, 0.5),
        "sharp_T": grid.base_raw,
        "random": projected_random,
        "random_near_threshold": near_threshold,
    }
    direct = {}
    for name, (field, direction) in _knot_fixtures(grid).items():
        projected_fixtures[name] = field
        if name == "checkerboard_saddle":
            direct[name] = _direction_check(area.evaluate, field, direction)
        else:
            # The bound is C1 and can have distinct curvatures at a knot.
            # Verify its derivative independently from each smooth side.
            direct[f"{name}_increasing"] = _direction_check(area.evaluate, field, direction, at_bound=True)
            direct[f"{name}_decreasing"] = _direction_check(area.evaluate, field, -direction, at_bound=True)
        # Compare both sides of the nonsmooth exact marching-squares event.
        projected_fixtures[f"{name}_below"] = field - 1e-4 * direction
        projected_fixtures[f"{name}_above"] = field + 1e-4 * direction
    direction = rng.uniform(-1.0, 1.0, grid.size)
    direct["random"] = _direction_check(area.evaluate, projected_random, direction)
    direct["random_near_threshold"] = _direction_check(area.evaluate, near_threshold, direction)
    comparisons = {name: _geometry_check(grid, area, field) for name, field in projected_fixtures.items()}
    empty = comparisons["all_solid"]["bound_mm2"]
    full = comparisons["all_fluid"]["bound_mm2"]
    all_ties = comparisons["all_threshold_ties"]["bound_mm2"]
    envelope_area = FIXED_FLUID_AREA_MM2 + grid.size * grid.cell_area_mm2
    uniform_errors = {"all_solid": abs(empty - FIXED_FLUID_AREA_MM2),
                      "all_fluid": abs(full - envelope_area),
                      "all_threshold_ties": abs(all_ties - envelope_area)}
    return {
        "grid_shape": grid.shape,
        "projected_field_checks": direct,
        "reconstruction_comparisons": comparisons,
        "uniform_envelope_area_errors_mm2": uniform_errors,
        "passed": all(item["passed"] for item in (*direct.values(), *comparisons.values()))
                  and max(uniform_errors.values()) <= EQUIVALENCE_TOLERANCE,
    }


def _archived_check(grid, archive, rng):
    checkpoint = json.loads(archive.read_text(encoding="utf-8"))
    if not checkpoint["completed"]:
        raise ValueError(f"Archived evaluation was incomplete: {archive}")
    saved = checkpoint["results"]
    raw = np.asarray(saved["design"], dtype=float)
    projected = np.asarray(saved["projected_density"], dtype=float)
    if raw.shape != (grid.size,) or projected.shape != raw.shape:
        raise ValueError("Archived design and projected field must fill the original design grid")
    stage = saved.get("stage", dict(q=0.01, beta=0.0, alpha_max=1e5))
    pipeline = DensityPipeline(grid)
    area = ReconstructedAreaUpperBound(grid, transition=TRANSITION, regularization=REGULARIZATION)

    def pulled_back(candidate, active_stage):
        density = pipeline.evaluate(candidate, **active_stage)
        value, covector = area.evaluate(density.projected)
        gradient = pipeline.filter_pullback(density.projection_derivative * covector)
        return value, gradient

    replayed = pipeline.evaluate(raw, **stage).projected
    reproduced_projected_error = float(np.max(np.abs(projected - replayed)))
    archived_geometry = _geometry_check(grid, area, projected)
    saved_area = saved.get("geometry", {}).get("physical_fluid_area_mm2")
    reproduced_geometry_error = (None if saved_area is None
                                 else abs(archived_geometry["reconstructed_mm2"] - saved_area))
    archived_geometry["archived_reconstructed_mm2"] = saved_area
    archived_geometry["archive_reconstructed_area_difference_mm2"] = reproduced_geometry_error
    direction = rng.uniform(0.25, 0.75, grid.size) * (1.0 - 2.0 * raw)
    raw_interior = 0.025 + 0.95 * raw
    beta4_stage = {**stage, "beta": 4.0}
    interior_direction = rng.uniform(-0.75, 0.75, grid.size)
    comparisons = {"archived_field": archived_geometry}
    for offset in (0.001, 0.005, 0.01):
        candidate = pipeline.evaluate(raw + offset * direction, **stage)
        comparisons[f"archive_inward_perturbation_{offset:.3f}"] = _geometry_check(
            grid, area, candidate.projected)
    projected_interior = pipeline.evaluate(raw_interior, **beta4_stage).projected
    comparisons["archive_interior_beta4"] = _geometry_check(grid, area, projected_interior)
    checks = {
        "archived_raw_direction": _direction_check(
            lambda candidate: pulled_back(candidate, stage), raw, direction, at_bound=True),
        "archived_nonlinear_projection_direction": _direction_check(
            lambda candidate: pulled_back(candidate, beta4_stage), raw_interior, interior_direction),
    }
    return {
        "archive_path": str(archive),
        "archive_stage": stage,
        "reproduced_projected_density_inf_error": reproduced_projected_error,
        "reconstruction_comparisons": comparisons,
        "raw_pullback_checks": checks,
        "passed": (reproduced_projected_error <= EQUIVALENCE_TOLERANCE
                   and (reproduced_geometry_error is None or
                        reproduced_geometry_error <= EQUIVALENCE_TOLERANCE)
                   and all(check["passed"] for check in (*checks.values(), *comparisons.values()))),
    }


def verify(archive):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run physical-area verification in a Sol Slurm job")
    rng = np.random.default_rng(43274)
    original = DesignGrid(0.5)
    cases = {"original_envelope": lambda: _envelope_check(original, rng),
             "expanded_envelope": lambda: _envelope_check(ExpandedDesignGrid(0.5), rng),
             "archived_model_b": lambda: _archived_check(original, archive, rng)}
    report = {}
    for name, check in cases.items():
        try:
            report[name] = check()
        except Exception:
            report[name] = {"passed": False, "error": traceback.format_exc()}
        print(f"{name}: {'pass' if report[name]['passed'] else 'fail'}", flush=True)
    return {
        "passed": all(case["passed"] for case in report.values()),
        "sol_job_id": os.environ["SLURM_JOB_ID"],
        "transition": TRANSITION,
        "regularization": REGULARIZATION,
        "limits": {"area_equivalence_mm2": EQUIVALENCE_TOLERANCE,
                   "scaled_derivative_error": DERIVATIVE_TOLERANCE,
                   "minimum_taylor_order": MINIMUM_TAYLOR_ORDER},
        "cases": report,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="New report directory")
    parser.add_argument("--archive", type=Path,
                        default=Path("/scratch") / os.environ.get("USER", "") / ARCHIVE_RELATIVE_PATH)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run physical-area verification in a Sol Slurm job")
    args.output.mkdir(parents=True, exist_ok=False)
    result = verify(args.archive)
    completed = write_report(args.output, "physical_area_verification", lambda: result)
    raise SystemExit(0 if completed and result["passed"] else 1)


if __name__ == "__main__":
    main()
