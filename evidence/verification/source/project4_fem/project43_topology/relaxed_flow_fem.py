"""Run the matched relaxed-split continuation on the original Taylor–Hood FEM."""

import argparse
from dataclasses import asdict, replace
from hashlib import sha256
import json
import os
from pathlib import Path

import numpy as np

from ..fem import State
from .connectivity import ConnectivityDiagnosticEvaluator
from .evaluator import REYNOLDS, Stage
from .full_density_continuation import MAX_AREA_MM2, ROUTES, restore_stage
from .gate2 import environment, write_report
from .gate2_restoration import FEASIBILITY_TOLERANCE, measured
from .optimize import optimize
from .reconstruction import reconstruct
from .relaxed_flow import (
    BANDS, MMA_BUDGET, RESTORE_BUDGET, BandedEvaluator, candidate_designs, independent_replay, prepare, save, setup,
)
from .restoration import restore


RECOVERY_AREA_BUFFERS_MM2 = (.25, .5, 1.)


def read_report(path):
    report = json.loads(path.read_text(encoding="utf-8"))
    if not report["completed"]:
        raise ValueError(f"Archive is incomplete: {path}")
    return report["results"]


def screen_history(hydraulic, records, output, name):
    """Reconstruct each feasible optimizer trial without running a PDE solve."""
    screened = []
    for index, row in enumerate(records):
        design = np.asarray(row["design"], dtype=float)
        if not (row["feasible"] and design.shape == (hydraulic.pipeline.grid.size,)
                and np.isfinite(design).all() and np.all((0. <= design) & (design <= 1.))):
            continue
        record = {"history_index": index, "objective": row["objective"], "design": design,
                  "inequalities": row["inequalities"], "exact_geometry_passed": False}
        try:
            geometry = reconstruct(hydraulic.pipeline.grid, hydraulic.density(design).projected)
        except (RuntimeError, ValueError) as error:
            record["error"] = f"{type(error).__name__}: {error}"
        else:
            record.update(area_mm2=geometry.area_mm2, connected_ports=geometry.port_component is not None,
                          exact_geometry_passed=geometry.port_component is not None
                          and geometry.area_mm2 <= MAX_AREA_MM2)
        screened.append(record)
    save(output, name, screened)
    return sorted((record for record in screened if record["exact_geometry_passed"]),
                  key=lambda record: record["objective"])


def replayable_trials(evaluator, candidates, maximum=4):
    """Select diverse contour-feasible starts with an independent original-constraint replay."""
    if not candidates:
        return []
    smallest = min(candidates, key=lambda row: row["area_mm2"])
    margin = min(candidates, key=lambda row: max(row["inequalities"][:2]))
    options = [candidates[0], smallest, margin, *candidates]
    accepted = []
    tried = set()
    for row in options:
        design = row["design"]
        key = design.tobytes()
        if key in tried:
            continue
        tried.add(key)
        try:
            evaluation = evaluator.evaluate(design)
            certificate, _ = evaluator.certificate(design)
        except (RuntimeError, ValueError, np.linalg.LinAlgError):
            continue
        if max(evaluation.constraints) <= FEASIBILITY_TOLERANCE and all(certificate.hard_paths):
            accepted.append((f"exact_history_{row['history_index']:03d}", design))
            if len(accepted) == maximum:
                break
    return accepted


def buffered_seed(evaluator, hydraulic, incoming, output):
    """Temporarily harden the projected-area cap during Phase I seed generation."""
    for index, buffer in enumerate(RECOVERY_AREA_BUFFERS_MM2):
        def constrained(design):
            evaluation = evaluator.evaluate(design)
            constraints = evaluation.constraints.copy()
            constraints[2] += buffer
            return replace(evaluation, constraints=constraints)

        restoration, cost = measured(lambda: restore(
            constrained, incoming, margin=0., max_evaluations=RESTORE_BUDGET, hard_tail=3), hydraulic)
        save(output, f"buffered_restoration_{index:02d}", {"buffer_mm2": buffer, "result": asdict(restoration), "cost": cost})
        history = [{"design": row.design, "objective": row.physical_objective,
                    "inequalities": row.inequalities, "feasible": row.engineering_feasible}
                   for row in restoration.trial_history]
        contour_candidates = screen_history(hydraulic, history, output, f"buffered_screen_{index:02d}")
        candidates = replayable_trials(evaluator, contour_candidates, maximum=1)
        if candidates:
            return candidates[0][1]
    return None


def run(output, model, common, recover_from=None):
    domain, pipeline, graph = setup()
    grid = pipeline.grid
    if grid.size != 240:
        raise RuntimeError("The original full-density comparison requires 240 design cells")
    domain.fem.inertia = REYNOLDS
    band = BANDS[model]
    common_reference = common / "baseline.npz"
    with np.load(common_reference) as archive:
        reference = float(archive["power"])
        values = archive["state"].copy()
    if not np.isfinite(reference) or reference <= 0. or values.shape != (domain.fem.size,) or not np.isfinite(values).all():
        raise ValueError("The common immersed reference does not match the original FEM")
    free = np.setdiff1d(np.arange(domain.fem.size), domain.fixed)
    seed = State(values, [], None, free)
    common_fitted = common / "sharp_fitted_baseline.json"
    fitted_archive = json.loads(common_fitted.read_text(encoding="utf-8"))
    if not fitted_archive["completed"] or not fitted_archive["results"]["fitted_diagnostic_passed"]:
        raise ValueError("The common sharp fitted reference failed its mesh or state check")
    reference_fitted_power = fitted_archive["results"]["levels"][-1]["mechanical_power_W_per_m"]
    if not np.isfinite(reference_fitted_power) or reference_fitted_power <= 0.:
        raise ValueError("The common fitted reference power must be positive")

    incoming = grid.base_raw.copy()
    stages = [Stage(.01, beta, alpha) for beta, alpha in ROUTES["coupled"]]
    settings = {
        "model": model, "flow_band": band, "inertia": REYNOLDS, "operator": "original Taylor–Hood FEM",
        "design_cells": grid.size, "design_spacing_mm": grid.spacing_mm, "fem_spacing_mm": .5,
        "initial_raw_sharp_t_design": incoming, "area_cap_mm2": MAX_AREA_MM2,
        "mma_budget": MMA_BUDGET, "restoration_budget": RESTORE_BUDGET,
        "stages": [asdict(stage) for stage in stages], "normalization_W_per_m": reference,
        "fitted_sharp_power_W_per_m": reference_fitted_power,
        "graph_cutoff": graph.cutoff, "graph_floor_budget": graph.floor_budget,
        "graph_minimum_conductance": graph.minimum_conductance,
        "common_reference_directory": str(common.resolve()),
        "common_immersed_sha256": sha256(common_reference.read_bytes()).hexdigest(),
        "common_fitted_sha256": sha256(common_fitted.read_bytes()).hexdigest(),
    }

    best = None
    completed = []
    if recover_from is not None:
        recover_from = recover_from.resolve()
        previous = read_report(recover_from / "summary.json")
        original = read_report(recover_from / "setup.json")
        keys = ("model", "flow_band", "common_immersed_sha256", "common_fitted_sha256", "stages")
        if any(original.get(key) != (list(settings[key]) if key == "flow_band" else settings[key]) for key in keys):
            raise ValueError("Recovery requires the same original lane, stage schedule, and common references")
        completed = previous["completed_stage_selections"].copy()
        if len(completed) >= len(stages):
            raise ValueError("The original lane already selected every continuation stage")
        best = previous["best_physical"]
        for index in range(len(completed)):
            stage_name = f"stage_{index:02d}"
            (output / stage_name).symlink_to(recover_from / stage_name, target_is_directory=True)
        if completed:
            checkpoint = read_report(recover_from / completed[-1] / "checkpoint.json")
            incoming = np.asarray(checkpoint["design"], dtype=float)
            previous_stage = stages[len(completed) - 1]
            evaluator = BandedEvaluator(domain, pipeline, previous_stage, reference, band=band, initial_state=seed)
            evaluator.evaluate(incoming)
            seed = evaluator.state
        settings.update(recovered_from=str(recover_from), recovery_stage_index=len(completed),
                        buffered_restoration_area_mm2=RECOVERY_AREA_BUFFERS_MM2)
    start_index = len(completed)
    save(output, "setup", settings)

    stop_reason = "completed schedule"
    for index, stage in enumerate(stages[start_index:], start_index):
        stage_dir = output / f"stage_{index:02d}"
        stage_dir.mkdir()
        hydraulic = BandedEvaluator(domain, pipeline, stage, reference, band=band, initial_state=seed)
        evaluator = ConnectivityDiagnosticEvaluator(hydraulic, graph)
        recovering = recover_from is not None and index == start_index
        candidates = []
        start = None
        if recovering:
            archive_dir = recover_from / stage_dir.name
            if (archive_dir / "restoration.json").is_file():
                restoration = read_report(archive_dir / "restoration.json")
                if restoration["accepted_design"] is not None:
                    start = np.asarray(restoration["accepted_design"], dtype=float)
            old_mma = archive_dir / "mma.json"
            if old_mma.is_file():
                archive = read_report(old_mma)["result"]
                ranked = screen_history(hydraulic, archive["evaluation_history"], stage_dir, "archived_exact_screen")
                candidates = replayable_trials(evaluator, ranked)
            if not candidates:
                start = buffered_seed(evaluator, hydraulic, incoming if start is None else start, stage_dir)
                save(stage_dir, "recovery_seed", {"stage": asdict(stage), "accepted_design": start})
        else:
            start, initial, restoration = restore_stage(evaluator, hydraulic, incoming, RESTORE_BUDGET, True)
            save(stage_dir, "restoration", {"stage": asdict(stage), "accepted_design": start,
                                           "initial": initial, "attempts": restoration})
        if candidates and start is None:
            start = incoming
        if start is None:
            stop_reason = "no feasible immersed start after restoration"
            break
        if not candidates:
            outcome, work = measured(lambda: optimize(evaluator, start, algorithm="mma",
                                                     max_evaluations=MMA_BUDGET), hydraulic)
            save(stage_dir, "mma", {"stage": asdict(stage), "result": asdict(outcome), "cost": work})
            candidates = candidate_designs(outcome, start)
            if recovering:
                rows = [asdict(point) for point in outcome.evaluation_history]
                ranked = screen_history(hydraulic, rows, stage_dir, "recovery_mma_exact_screen")
                candidates = replayable_trials(evaluator, ranked) + candidates
        else:
            save(stage_dir, "selected_archived_trials", {"source": str(archive_dir), "designs": candidates})
        audited = []
        audited_designs = set()
        for label, design in candidates:
            key = np.asarray(design).tobytes()
            if key in audited_designs:
                continue
            audited_designs.add(key)
            candidate_dir = stage_dir / label
            candidate_dir.mkdir()
            replay = BandedEvaluator(domain, pipeline, stage, reference, band=band, initial_state=seed)
            try:
                report, audit_work = measured(
                    lambda: independent_replay(replay, graph, design, candidate_dir, reference_fitted_power), replay)
                density = replay.density(design)
            except (RuntimeError, ValueError, AssertionError, np.linalg.LinAlgError) as error:
                save(candidate_dir, "checkpoint", {"design": design, "stage": asdict(stage),
                                                    "error": f"{type(error).__name__}: {error}"})
                continue
            report.update(stage=asdict(stage), projected_density=density.projected,
                          filtered_density=density.filtered, audit_cost=audit_work)
            save(candidate_dir, "checkpoint", report)
            relative = candidate_dir.relative_to(output).as_posix()
            if report["fitted_passed"]:
                qualified = {
                    "checkpoint": relative, "stage": asdict(stage),
                    "fitted_power_W_per_m": report["fitted_power_W_per_m"],
                    "fitted_power_ratio_vs_sharp": report["fitted_power_ratio_vs_sharp"],
                    "fitted_side_share": report["fitted"]["levels"][-1]["side_share"],
                    "fitted_area_mm2": report["geometry"]["physical_fluid_area_mm2"],
                    "optimizer_sample": not np.array_equal(design, start),
                    "kkt_passed": report["kkt"]["verified"],
                }
                if best is None or (not qualified["optimizer_sample"], qualified["fitted_power_W_per_m"]) < (
                        not best["optimizer_sample"], best["fitted_power_W_per_m"]):
                    best = qualified
                save(candidate_dir, "physical_acceptance", qualified)
            if (report["immersed_feasible"] and report["geometry"]["connected_ports"]
                    and report["geometry"]["physical_fluid_area_mm2"] <= MAX_AREA_MM2):
                audited.append((report, replay.state, relative))
        if not audited:
            stop_reason = "no independently replayed connected and area-feasible handoff"
            break
        report, seed, path = min(audited, key=lambda row: (not row[0]["fitted_passed"], row[0]["objective"]))
        incoming = report["design"]
        completed.append(path)
        save(stage_dir, "selected", {"stage": asdict(stage), "checkpoint": path, "best_physical_so_far": best})
        print(f"{model} original FEM stage {index} complete, best physical {best}", flush=True)
    return {"model": model, "flow_band": band, "completed_stage_selections": completed,
            "scheduled_stages": len(stages), "best_physical": best, "stop_reason": stop_reason,
            "recovered_from": None if recover_from is None else str(recover_from)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", *BANDS))
    parser.add_argument("output", type=Path)
    parser.add_argument("--common", type=Path)
    parser.add_argument("--recover-from", type=Path, help="Completed original FEM lane stopped before selecting all stages")
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical study as a Sol Slurm job")
    if args.mode != "prepare" and args.common is None:
        parser.error("Supply the completed common baseline directory")
    if args.mode == "prepare" and args.recover_from is not None:
        parser.error("Recovery applies to a model lane")
    args.output.mkdir(parents=True, exist_ok=False)
    write_report(args.output, "environment", environment)
    if not write_report(args.output, "summary", lambda: (
            prepare(args.output) if args.mode == "prepare" else run(args.output, args.mode, args.common, args.recover_from))):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
