"""Recover one common hybrid diagnostic seed from an archived baffle sweep."""

import argparse
from dataclasses import asdict, replace
import json
import os
from pathlib import Path

import numpy as np

from .connected_split_diagnostic import SplitRestoration
from .connectivity import ConnectivityGraph, RectangularBaffleClearance
from .density import DesignGrid
from .domain import build_domain
from .gate2 import environment, write_report
from .gate2_restoration import FEASIBILITY_TOLERANCE, measured
from .hybrid_density import BAFFLE_X_MM, HEIGHT_LIMIT, HEIGHT_MM, HybridDensityModel, HybridEvaluator, hybrid_reference
from .hybrid_driver import (
    BETA_ZERO_STAGE, MAX_AREA_MM2, MIN_SEED_HEIGHT_MM, UnitBoxEvaluator,
    audit_candidate, exact_feasible, physical_design, preserve, unit_design,
)
from .optimize import optimize


ALPHA_MAX_CHOICES = (1e5, 3e5, 1e6)


def read_report(path):
    report = json.loads(path.read_text(encoding="utf-8"))
    if not report["completed"]:
        raise ValueError(f"Incomplete archived report: {path}")
    return report["results"]


def archived_starts(sweep, model, maximum):
    eligible = []
    for row in sweep:
        height = row["height_mm"]
        immersed = row.get("hybrid_immersed")
        fitted = row.get("fitted_sharp_baffle", {})
        paths = row.get("strict_graph_paths", ())
        if (not MIN_SEED_HEIGHT_MM < height <= HEIGHT_MM * HEIGHT_LIMIT
                or not row.get("reconstructed_connected")
                or len(paths) != 2 or not all(paths)
                or row["reconstructed_area_mm2"] > MAX_AREA_MM2
                or immersed is None or max(immersed["inequalities"][2:]) > FEASIBILITY_TOLERANCE
                or not fitted.get("state_passed")):
            continue
        expected = model.initial.copy()
        expected[-1] = height / HEIGHT_MM
        if not np.array_equal(np.asarray(row["physical_design"]), expected):
            raise ValueError("Archived baffle sweep does not match the original sharp-T hybrid preset")
        eligible.append(row)

    immersed_order = sorted(eligible, key=lambda row: abs(row["hybrid_immersed"]["side_fraction"] - .5))
    fitted_order = sorted(eligible, key=lambda row: abs(row["fitted_sharp_baffle"]["side_share"] - .5))
    order = immersed_order[:1] + fitted_order[:1] + immersed_order
    distinct = {}
    for row in order:
        distinct.setdefault(row["height_mm"], row)
    return list(distinct.values())[:maximum]


def search(output, prior, algorithm, budget, maximum, alpha_max):
    sweep = read_report(prior / "baffle_only_sweep.json")
    setup = read_report(prior / "setup.json")
    if not setup["stages"] or setup["stages"][0] != asdict(BETA_ZERO_STAGE):
        raise ValueError("Archived baffle sweep has a different initial continuation stage")
    stage = replace(BETA_ZERO_STAGE, alpha_max=alpha_max)
    grid = DesignGrid(.5)
    model = HybridDensityModel(grid)
    starts = archived_starts(sweep, model, maximum)
    if len(starts) < 2:
        return {"successful": False, "shared_unit_design": None, "shared_physical_design": None,
                "attempted_initial_heights_mm": [], "stage": asdict(stage),
                "reason": "fewer than two distinct eligible archived heights"}
    domain = build_domain(grid, .5)
    graph = ConnectivityGraph(grid, baffle=RectangularBaffleClearance(
        x_min_mm=BAFFLE_X_MM[0], x_max_mm=BAFFLE_X_MM[1], max_height_mm=HEIGHT_MM * HEIGHT_LIMIT))
    reference, state = hybrid_reference(domain, model, BETA_ZERO_STAGE.alpha_max)
    np.testing.assert_allclose(reference, setup["normalized_power_reference_W_per_m"], rtol=1e-10, atol=0.)
    preserve(output, "setup", {"source_diagnostic_directory": str(prior.resolve()),
                               "source_seed_job_id": setup["slurm_job_id"],
                               "stage": asdict(stage), "algorithm": algorithm,
                               "power_reference_alpha_max": BETA_ZERO_STAGE.alpha_max,
                               "normalized_power_reference_W_per_m": reference,
                               "evaluations_per_start": budget, "ranked_starts": starts,
                               "start_ranking_alpha_max": BETA_ZERO_STAGE.alpha_max})

    def fresh_evaluator():
        return HybridEvaluator(domain, model, stage, reference, initial_state=state, graph=graph)

    attempts, selected = [], None
    for index, row in enumerate(starts):
        folder = output / f"start_{index:02d}"
        folder.mkdir()
        hybrid = fresh_evaluator()
        wrapper = SplitRestoration(UnitBoxEvaluator(hybrid))
        initial = unit_design(row["physical_design"])
        outcome, cost = measured(lambda: optimize(
            wrapper, initial, algorithm=algorithm, max_evaluations=budget), hybrid)
        preserve(folder, "split_search", {"result": asdict(outcome), "cost": cost})

        feasible = sorted(outcome.feasible_history, key=lambda point: point.objective)
        checks = []
        for point in feasible:
            # The split search objective is the squared distance to 50 percent.
            if point.objective > (.001 * (1. + FEASIBILITY_TOLERANCE))**2 + 1e-15:
                break
            try:
                accepted = exact_feasible(fresh_evaluator(), point.design, positive=True)
                check = {"unit_design": point.design, "split_objective": point.objective, "accepted": accepted}
            except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
                check = {"unit_design": point.design, "split_objective": point.objective,
                         "accepted": False, "error": f"{type(error).__name__}: {error}"}
            checks.append(check)
            if check["accepted"]:
                selected = point.design.copy()
                break

        audited = selected if selected is not None else outcome.best_feasible_design
        if audited is not None:
            audit_dir = folder / "replayed_candidate"
            audit_dir.mkdir()
            write_report(audit_dir, "checkpoint", lambda: audit_candidate(fresh_evaluator(), audited, audit_dir)[0])
        preserve(folder, "exact_replay", checks)
        attempts.append({"height_mm": row["height_mm"], "directory": folder.name,
                         "best_area_and_graph_objective": feasible[0].objective if feasible else None,
                         "replayed_split_band_trials": len(checks)})
        if selected is not None:
            break

    return {"successful": selected is not None, "shared_unit_design": selected,
            "shared_physical_design": None if selected is None else physical_design(selected),
            "attempted_initial_heights_mm": [attempt["height_mm"] for attempt in attempts],
            "stage": asdict(stage), "attempts": attempts, "production_verified": False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prior_diagnostic", type=Path, help="First seed job's diagnostic directory")
    parser.add_argument("output", type=Path, help="New diagnostic directory for the shared seed")
    parser.add_argument("--algorithm", choices=("mma", "slsqp"), default="mma")
    parser.add_argument("--alpha-max", type=float, choices=ALPHA_MAX_CHOICES, default=BETA_ZERO_STAGE.alpha_max)
    parser.add_argument("--max-evaluations", type=int, default=150)
    parser.add_argument("--starts", type=int, choices=(2, 3), default=3)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical diagnostic as a Sol Slurm job")
    if args.max_evaluations < 1:
        parser.error("The evaluation budget must be positive")
    args.output.mkdir(parents=True, exist_ok=False)
    preserve(args.output, "environment", environment())
    if not write_report(args.output, "seed", lambda: search(
            args.output, args.prior_diagnostic, args.algorithm, args.max_evaluations, args.starts, args.alpha_max)):
        raise SystemExit(1)
    if not read_report(args.output / "seed.json")["successful"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
