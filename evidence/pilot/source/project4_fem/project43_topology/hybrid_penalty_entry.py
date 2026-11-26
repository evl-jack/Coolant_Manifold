"""Run matched hybrid model diagnostics from a verified high-penalty shared seed."""

import argparse
from dataclasses import asdict
from hashlib import sha256
import json
import os
from pathlib import Path

import numpy as np

from .connectivity import ConnectivityGraph, RectangularBaffleClearance
from .density import DesignGrid
from .domain import build_domain
from .evaluator import Stage
from .gate2 import environment, write_report
from .hybrid_density import BAFFLE_X_MM, HEIGHT_LIMIT, HEIGHT_MM, HybridDensityModel, HybridEvaluator, hybrid_reference
from .hybrid_driver import (
    BETA_ZERO_STAGE, HEIGHT_COSTS, ROUTES, audit_candidate, exact_feasible, physical_design, preserve, run_model,
)


ENTRY_SCHEDULES = {
    3e5: ((0., 3e5), (1., 3e5), (2., 3e5), (2., 1e6), (4., 1e6), (8., 1e6)),
    1e6: ((0., 1e6), (1., 1e6), (2., 1e6), (4., 1e6), (8., 1e6)),
}


def run(output, seed_path, mode, max_stages, mma_budget, restore_budget):
    seed_bytes = seed_path.read_bytes()
    archived = json.loads(seed_bytes)
    if not archived["completed"] or not archived["results"]["successful"]:
        raise ValueError("The archived diagnostic did not establish a common seed")
    seed = archived["results"]
    seed_stage = Stage(**seed["stage"])
    if seed_stage.q != BETA_ZERO_STAGE.q or seed_stage.beta != 0. or seed_stage.alpha_max not in ENTRY_SCHEDULES:
        raise ValueError("The seed must use q=0.01, beta=0 and alpha_max in {3e5, 1e6}")
    stages = tuple(Stage(seed_stage.q, beta, alpha) for beta, alpha in ENTRY_SCHEDULES[seed_stage.alpha_max])
    if stages[0] != seed_stage:
        raise ValueError("The first scheduled stage must reproduce the archived seed stage")
    lower_penalty_stages = [asdict(Stage(seed_stage.q, beta, alpha))
                            for beta, alpha in ROUTES["coupled"] if alpha < seed_stage.alpha_max]
    grid = DesignGrid(.5)
    model = HybridDensityModel(grid)
    common = np.asarray(seed["shared_unit_design"], dtype=float)
    if common.shape != (model.size,) or not np.isfinite(common).all() or np.any((common < 0.) | (common > 1.)):
        raise ValueError("The archived seed has invalid unit-box coordinates")
    np.testing.assert_array_equal(physical_design(common), seed["shared_physical_design"])
    domain = build_domain(grid, .5)
    graph = ConnectivityGraph(grid, baffle=RectangularBaffleClearance(
        x_min_mm=BAFFLE_X_MM[0], x_max_mm=BAFFLE_X_MM[1], max_height_mm=HEIGHT_MM * HEIGHT_LIMIT))
    reference, sharp_state = hybrid_reference(domain, model, BETA_ZERO_STAGE.alpha_max)
    scheduled = stages[:max_stages]
    setup = {
        "diagnostic_only": True, "slurm_job_id": os.environ["SLURM_JOB_ID"],
        "shared_seed_path": str(seed_path.resolve()), "shared_seed_sha256": sha256(seed_bytes).hexdigest(),
        "shared_unit_design": common, "shared_physical_design": physical_design(common), "seed_stage": asdict(seed_stage),
        "legacy_sharp_reference_alpha_max": BETA_ZERO_STAGE.alpha_max, "normalized_power_reference_W_per_m": reference,
        "design_spacing_mm": grid.spacing_mm, "design_variables": model.size, "mesh_triangles": len(domain.fem.triangles),
        "planned_entry_stages": [asdict(stage) for stage in stages], "scheduled_entry_stages": [asdict(stage) for stage in scheduled],
        "lower_penalty_stages_skipped": lower_penalty_stages,
        "reason_for_lower_penalty_skip": "Independent high-penalty entry diagnostic at the archived shared-seed stage",
        "height_costs_by_model": HEIGHT_COSTS, "mma_budget_per_stage": mma_budget,
        "restoration_budget_per_attempt": restore_budget,
    }
    preserve(output, "setup", setup)

    results = {}
    for model_name in (tuple(HEIGHT_COSTS) if mode == "all" else (mode,)):
        lane = output / model_name
        lane.mkdir()
        replay = lane / "seed_replay"
        replay.mkdir()

        def run_lane():
            evaluator = HybridEvaluator(domain, model, seed_stage, reference, initial_state=sharp_state,
                                        height_cost=HEIGHT_COSTS[model_name], graph=graph)
            feasible = exact_feasible(evaluator, common, positive=True)
            preserve(replay, "acceptance", {"stage": asdict(seed_stage), "unit_design": common,
                                            "positive_exact_feasible": feasible})
            if not feasible:
                raise ValueError("The common seed failed the independent replay at its archived penalty")
            audit, _ = audit_candidate(evaluator, common, replay)
            preserve(replay, "checkpoint", audit)
            if not audit["exact_and_immersed_feasible"]:
                raise ValueError("The common seed failed independently integrated flow or exact reconstruction")
            context = domain, model, reference, evaluator.state, graph
            return run_model(lane, context, model_name, common.copy(), scheduled, mma_budget, restore_budget)

        completed = write_report(lane, "summary", run_lane)
        results[model_name] = {"completed": completed, "report": f"{model_name}/summary.json"}
    return {"diagnostic_only": True, "common_seed_stage": asdict(seed_stage), "model_results": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shared_seed", type=Path, help="Accepted high-penalty seed.json")
    parser.add_argument("output", type=Path, help="New per-job evidence directory")
    parser.add_argument("--mode", choices=("model-a", "model-b", "all"), required=True)
    parser.add_argument("--max-stages", type=int)
    parser.add_argument("--mma-evaluations", type=int, default=150)
    parser.add_argument("--restoration-evaluations", type=int, default=250)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical diagnostic as a Sol Slurm job")
    if min(args.mma_evaluations, args.restoration_evaluations) < 1 or (args.max_stages is not None and args.max_stages < 1):
        parser.error("Stage and evaluation budgets must be positive")
    args.output.mkdir(parents=True, exist_ok=False)
    preserve(args.output, "environment", environment())
    if not write_report(args.output, "summary", lambda: run(
            args.output, args.shared_seed, args.mode, args.max_stages, args.mma_evaluations, args.restoration_evaluations)):
        raise SystemExit(1)
    report = json.loads((args.output / "summary.json").read_text(encoding="utf-8"))["results"]
    if any(not lane["completed"] for lane in report["model_results"].values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
