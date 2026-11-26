"""Run the archived fixed-pipe continuation with the linear Stokes state."""

import argparse
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
from time import perf_counter

import numpy as np

from ..evaluator import REYNOLDS, Stage
from ..fixed_pipe_pilot import choose_checkpoint, progress_callback, save, snapshot
from ..fixed_pipe_pilot_model import SIDE_CONSTRAINT_SCALE, original_geometry, reconstruct_geometry
from ..gate2 import environment, write_report
from ..optimize import optimize
from ..restoration import FEASIBILITY_TOLERANCE, restore
from ..stokes_fixed_pipe_model import StokesFixedPipeEvaluator, build_stokes_fixed_pipe
from .fixed_pipe_kkt import TOLERANCES, archived_result, audit_point


def stage_progress(journal, history, evaluator, stage, phase, started):
    write_progress = progress_callback(journal, history, evaluator, stage, phase, started)
    last_unique = 0

    def progress(point, unique, total):
        nonlocal last_unique
        write_progress(point, unique, total)
        if unique == last_unique:
            return
        last_unique = unique
        design = np.asarray(point.design)
        if evaluator.last_design is None or not np.array_equal(design, evaluator.last_design):
            raise RuntimeError("The optimizer callback does not match the evaluated Stokes design")
        audited = audit_point(design, evaluator.last_evaluation, evaluator)
        row = {
            "stage": stage, "phase": phase, "unique_evaluations": unique,
            **{key: audited[key] for key in ("stationarity_inf", "primal_violation", "dual_violation",
                                             "complementarity_inf", "verified")},
            "multipliers": audited["multipliers"].tolist(),
        }
        print(json.dumps(row, allow_nan=False), file=journal, flush=True)

    return progress


def run_stage(folder, source, domain, pipeline, stage, number, incoming, budget, continuation_budget, journal, history):
    started = perf_counter()
    folder.mkdir()
    pilot_metrics_path = source / "original" / "metrics.json"
    pilot_metrics = archived_result(pilot_metrics_path)
    power_reference = pilot_metrics["power_w_per_m"]
    side_reference = pilot_metrics["side_fraction"]
    reference = StokesFixedPipeEvaluator(domain, pipeline, stage, power_reference, side_reference)
    original = pipeline.initial
    baseline = reference.evaluate(original)
    baseline_snapshot = snapshot(folder / "original", reference, original, baseline, original_geometry())
    evaluator = StokesFixedPipeEvaluator(domain, pipeline, stage, power_reference, side_reference,
                                         initial_state=reference.state)
    at_start = evaluator.evaluate(incoming)
    initial_g = float(at_start.constraints[0])
    print(f"Stage {number}, beta={stage.beta:g}: Stokes original side={baseline.side_fraction:.8f}, "
          f"required side={side_reference:.8f}, incoming constraint={initial_g:.6g}", flush=True)
    report = {
        "stage": number, "settings": asdict(stage), "physical_reynolds": REYNOLDS, "pde_inertia": 0.,
        "baseline_reference_model": "archived Navier–Stokes at this continuation stage",
        "ns_reference_sha256": hashlib.sha256(pilot_metrics_path.read_bytes()).hexdigest(),
        "power_reference_w_per_m": power_reference, "side_reference": side_reference,
        "incoming_side_fraction": at_start.side_fraction, "incoming_constraint": initial_g,
        "original": baseline_snapshot["metrics"], "original_satisfies_ns_side_floor":
        bool(baseline.constraints[0] <= FEASIBILITY_TOLERANCE),
        "reference_runtime_seconds": reference.forward_seconds + reference.adjoint_seconds,
        "restoration_attempted": initial_g > FEASIBILITY_TOLERANCE,
        "restoration_unique_evaluations": 0, "restarted_at_original": False,
    }
    history.append({"stage": number, "phase": "reference", "unique_evaluations": 0,
                    "normalized_power": baseline.objective, "side_fraction": baseline.side_fraction,
                    "reference_side_fraction": side_reference})
    start = incoming
    if report["restoration_attempted"]:
        restoration = restore(evaluator, incoming, margin=0., max_evaluations=budget,
                              progress=progress_callback(journal, history, evaluator, number, "restoration", started))
        save(folder, "restoration", asdict(restoration))
        report["restoration_unique_evaluations"] = restoration.unique_evaluations
        if restoration.accepted_design is not None:
            start = np.asarray(restoration.accepted_design)
        elif report["original_satisfies_ns_side_floor"]:
            start = original
            report["restarted_at_original"] = True
        else:
            report.update(feasible_checkpoint=False, elapsed_seconds=perf_counter() - started,
                          failure="Bounded restoration found no feasible design at the archived side-flow floor.")
            save(folder, "summary", report)
            return None, report
    result = optimize(evaluator, start, algorithm="MMA", max_evaluations=budget, mma_verbosity=0,
                      progress=stage_progress(journal, history, evaluator, number, "mma", started))
    save(folder, "mma", asdict(result))
    best_from, selected, value = choose_checkpoint(evaluator, result, start)
    initial_audit = audit_point(selected, value, evaluator)
    save(folder, "kkt_initial", initial_audit)
    extra_count = 0
    if continuation_budget and not initial_audit["verified"]:
        initial_geometry = reconstruct_geometry(pipeline, evaluator.density(selected))
        snapshot(folder / "initial_selected", evaluator, selected, value, initial_geometry)
        extra = optimize(evaluator, selected, algorithm="MMA", max_evaluations=continuation_budget, mma_verbosity=0,
                         progress=stage_progress(journal, history, evaluator, number, "mma_continuation", started))
        save(folder, "mma_continuation", asdict(extra))
        extra_count = extra.unique_evaluations
        best_from, selected, value = choose_checkpoint(evaluator, extra, selected)
        best_from = f"final-stage continuation: {best_from}"
    geometry = reconstruct_geometry(pipeline, evaluator.density(selected))
    selected_snapshot = snapshot(folder / "selected", evaluator, selected, value, geometry)
    kkt = audit_point(selected, value, evaluator)
    save(folder, "kkt", {**kkt, "settings": asdict(stage), "tolerances": TOLERANCES})
    report.update(
        feasible_checkpoint=True, selected=selected_snapshot["metrics"], selected_from=best_from,
        elapsed_seconds=perf_counter() - started, forward_seconds=evaluator.forward_seconds,
        adjoint_seconds=evaluator.adjoint_seconds, forward_solves=evaluator.forward_solves,
        adjoint_rhs=evaluator.adjoint_rhs, mma_native_status=result.status, mma_native_message=result.message,
        mma_unique_evaluations=result.unique_evaluations, mma_objective_calls=result.objective_calls,
        mma_evaluator_requests=result.evaluator_requests, extra_mma_unique_evaluations=extra_count,
        kkt_verified=kkt["verified"], selected_constraint=float(value.constraints[0]),
    )
    save(folder, "summary", report)
    print(f"Stage {number}: MMA {result.status} ({result.message}); power={value.power_w_per_m:.9g} W/m, "
          f"side={value.side_fraction:.8f}, KKT={kkt['verified']}", flush=True)
    return selected, report


def run(source, output, continuation_budget):
    started = perf_counter()
    archived_configuration = archived_result(source / "configuration.json")
    stages = [Stage(**parameters) for parameters in archived_configuration["stages"]]
    budget = archived_configuration["unique_evaluation_limit"]
    domain, pipeline = build_stokes_fixed_pipe()
    fem = domain.fem
    with np.load(source / "mesh.npz", allow_pickle=False) as mesh:
        for name, actual in (("points", fem.points), ("triangles", fem.triangles),
                             ("velocity_dofs", fem.velocity_dofs), ("cell_index", domain.cell_index)):
            np.testing.assert_array_equal(mesh[name], actual)
    original_checkpoint = source / "stage_01_beta_0" / "original" / "fields.npz"
    with np.load(original_checkpoint, allow_pickle=False) as archive:
        np.testing.assert_array_equal(archive["design"], pipeline.initial)
    save(output, "configuration", {
        **archived_configuration, "environment": environment(), "sol_job_id": os.environ.get("SLURM_JOB_ID"),
        "physical_reynolds": REYNOLDS, "pde_inertia": fem.inertia, "side_constraint_scale": SIDE_CONSTRAINT_SCALE,
        "power_reference": "Archived Navier–Stokes unchanged immersed original at the same continuation stage",
        "side_constraint": "Side fraction >= archived Navier–Stokes unchanged immersed original at the same stage",
        "final_stage_continuation_budget": continuation_budget,
    })
    np.savez_compressed(output / "mesh.npz", points=fem.points, triangles=fem.triangles,
                        velocity_points=fem.velocity_points, velocity_dofs=fem.velocity_dofs,
                        cell_index=domain.cell_index, density_x_mm=pipeline.grid.x_edges_mm,
                        density_y_mm=pipeline.grid.y_edges_mm)
    design = pipeline.initial
    history, reports = [], []
    with (output / "progress.jsonl").open("x", encoding="utf-8") as journal:
        for number, stage in enumerate(stages, start=1):
            name = f"stage_{number:02d}_beta_{stage.beta:g}"
            design, report = run_stage(output / name, source / name, domain, pipeline, stage, number, design,
                                       budget, continuation_budget if number == len(stages) else 0, journal, history)
            reports.append(report)
            if design is None:
                break
    return {"stages": reports, "elapsed_seconds": perf_counter() - started,
            "all_stages_feasible": len(reports) == len(stages) and all(r["feasible_checkpoint"] for r in reports),
            "final_kkt_verified": bool(reports[-1].get("kkt_verified")) if len(reports) == len(stages) else False,
            "mma_unique_evaluations_total": sum(r.get("mma_unique_evaluations", 0) for r in reports),
            "extra_mma_unique_evaluations": sum(r.get("extra_mma_unique_evaluations", 0) for r in reports),
            "forward_seconds": sum(r.get("forward_seconds", 0.) for r in reports),
            "adjoint_seconds": sum(r.get("adjoint_seconds", 0.) for r in reports)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archived_pilot", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--continuation-budget", type=int, default=100)
    args = parser.parse_args()
    if args.continuation_budget < 0:
        parser.error("--continuation-budget must be nonnegative")
    args.output.mkdir(parents=True, exist_ok=False)
    complete = write_report(args.output, "summary", lambda: run(args.archived_pilot, args.output, args.continuation_budget))
    status = archived_result(args.output / "summary.json") if complete else {}
    raise SystemExit(0 if status.get("all_stages_feasible") and status.get("final_kkt_verified") else 1)


if __name__ == "__main__":
    main()
