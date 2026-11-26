"""Single exploratory fixed-pipe Brinkman transition run."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
import json
import os
from pathlib import Path
from time import perf_counter

import numpy as np

from .evaluator import REYNOLDS, Stage
from .fixed_pipe_pilot_model import (
    SIDE_CONSTRAINT_SCALE, FixedPipeEvaluator, build_fixed_pipe, metrics_for, original_geometry, reconstruct_geometry,
)
from .fixed_pipe_pilot_plots import plot_stage
from .gate2 import environment, write_report
from .optimize import optimize
from .reconstruction import write_reconstruction
from .restoration import FEASIBILITY_TOLERANCE, restore


STAGES = tuple(Stage(.10, beta, 1e5) for beta in (0., 2., 8.))
BUDGET = 100


def save(output, name, value):
    if not write_report(output, name, lambda: value):
        raise RuntimeError(f"Could not save {name} to {output}")


def snapshot(folder, evaluator, design, value, geometry):
    folder.mkdir()
    density = evaluator.density(design)
    state = evaluator.state.values.copy()
    metrics = metrics_for(value, geometry, evaluator.pipeline)
    metrics["normalized_power"] = float(value.power_w_per_m / evaluator.power_reference)
    metrics["pumping_power_change_percent"] = 100. * (metrics["normalized_power"] - 1.)
    metrics["side_fraction"] = float(value.side_fraction)
    metrics["main_fraction"] = float(value.flux[1] / value.flux[0])
    write_reconstruction(geometry, folder / "geometry.json")
    np.savez_compressed(folder / "fields.npz", state=state, design=design, raw=density.raw,
                        filtered=density.filtered, projected=density.projected, resistance=density.resistance)
    save(folder, "metrics", metrics)
    return {"geometry": geometry, "density": density, "state": state, "metrics": metrics}


def progress_callback(journal, history, evaluator, stage, phase, started):
    last_unique = 0

    def progress(point, unique, total):
        nonlocal last_unique
        if unique == last_unique:
            return
        last_unique = unique
        normalized = point.objective if phase == "mma" else point.physical_objective
        side = evaluator.side_reference - point.inequalities[0] * SIDE_CONSTRAINT_SCALE
        row = {"stage": stage, "phase": phase, "unique_evaluations": unique, "total_requests": total,
               "elapsed_seconds": perf_counter() - started, "normalized_power": float(normalized),
               "side_fraction": float(side), "reference_side_fraction": float(evaluator.side_reference),
               "constraint": float(point.inequalities[0]), "forward_seconds": evaluator.forward_seconds,
               "adjoint_seconds": evaluator.adjoint_seconds}
        if phase == "restoration":
            row["epigraph_slack"] = point.slack
        history.append(row)
        line = json.dumps(row, allow_nan=False)
        print(line, file=journal, flush=True)
        print(line, flush=True)

    return progress


def choose_checkpoint(evaluator, result, start):
    """Replay the best feasible evaluated point, retaining native status separately."""
    candidates = (("best feasible MMA evaluation", result.best_feasible_design),
                  ("native terminal", result.terminal_design), ("stage start", start))
    for label, design in candidates:
        if design is None:
            continue
        try:
            value = evaluator.evaluate(design)
        except (RuntimeError, ValueError, np.linalg.LinAlgError):
            continue
        if value.constraints[0] <= FEASIBILITY_TOLERANCE:
            return label, design.copy(), value
    raise RuntimeError("No feasible stage checkpoint could be replayed")


def run_stage(folder, domain, pipeline, stage, stage_no, incoming, ref_state, journal, history):
    folder.mkdir()
    stage_start = perf_counter()
    original = pipeline.initial
    reference = FixedPipeEvaluator(domain, pipeline, stage, 1., 0., initial_state=ref_state)
    baseline_value = reference.evaluate(original)
    evaluator = FixedPipeEvaluator(domain, pipeline, stage, baseline_value.power_w_per_m,
                                   baseline_value.side_fraction, initial_state=reference.state)
    reference.power_reference = baseline_value.power_w_per_m
    baseline = snapshot(folder / "original", reference, original, baseline_value, original_geometry())
    starting_value = evaluator.evaluate(incoming)
    print(f"Stage {stage_no}, beta={stage.beta:g}: baseline side={baseline_value.side_fraction:.8f}, "
          f"continued start side={starting_value.side_fraction:.8f}", flush=True)
    history.append({"stage": stage_no, "phase": "reference", "unique_evaluations": 0,
                    "normalized_power": 1., "side_fraction": baseline_value.side_fraction,
                    "reference_side_fraction": baseline_value.side_fraction})
    restoration = None
    start = incoming
    if starting_value.constraints[0] > FEASIBILITY_TOLERANCE:
        restoration = restore(evaluator, incoming, margin=0., max_evaluations=BUDGET, progress=progress_callback(
            journal, history, evaluator, stage_no, "restoration", stage_start))
        save(folder, "restoration", asdict(restoration))
        # The unchanged design is a known feasible restart at every stage.
        start = original if restoration.accepted_design is None else np.asarray(restoration.accepted_design)
    result = optimize(evaluator, start, algorithm="MMA", max_evaluations=BUDGET, mma_verbosity=0, progress=progress_callback(
        journal, history, evaluator, stage_no, "mma", stage_start))
    save(folder, "mma", asdict(result))
    selected_from, design, value = choose_checkpoint(evaluator, result, start)
    geometry = reconstruct_geometry(pipeline, evaluator.density(design))
    candidate = snapshot(folder / "selected", evaluator, design, value, geometry)
    report = {"stage": stage_no, "settings": asdict(stage), "elapsed_seconds": perf_counter() - stage_start,
              "reference_runtime_seconds": reference.forward_seconds + reference.adjoint_seconds,
              "forward_seconds": evaluator.forward_seconds, "adjoint_seconds": evaluator.adjoint_seconds,
              "incoming_side_fraction": starting_value.side_fraction, "restoration_attempted": restoration is not None,
              "restoration_unique_evaluations": 0 if restoration is None else restoration.unique_evaluations,
              "restarted_at_original": restoration is not None and restoration.accepted_design is None,
              "mma_native_status": result.status, "mma_native_message": result.message,
              "mma_unique_evaluations": result.unique_evaluations, "mma_objective_calls": result.objective_calls,
              "mma_evaluator_requests": result.evaluator_requests, "selected_from": selected_from,
              "original": baseline["metrics"], "selected": candidate["metrics"]}
    save(folder, "summary", report)
    print(f"Stage {stage_no} selected {value.objective:.9f}; MMA {result.status}: {result.message}", flush=True)
    return design, reference.state, baseline, candidate, report


def run(output):
    started = perf_counter()
    domain, pipeline = build_fixed_pipe()
    design, ref_state = pipeline.initial, None
    fem = domain.fem
    job = os.environ.get("SLURM_JOB_ID")
    save(output, "configuration", {"sol_job_id": job, "environment": environment(),
        "stages": [asdict(stage) for stage in STAGES], "unique_evaluation_limit": BUDGET,
        "reynolds": REYNOLDS, "density_spacing_mm": .25, "fem_spacing_mm": domain.max_spacing_mm,
        "filter_radius_mm": .5, "reconstruction_threshold": .5, "design_rectangle_mm": [8., 4., 12., 6.],
        "main_pipe_mm": [0., 0., 20., 4.], "outlet_collar_mm": [9.5, 6., 10.5, 7.],
        "editable_density_variables": int(design.size), "triangles": len(fem.triangles),
        "power_reference": "Unchanged immersed original at the same continuation stage",
        "side_constraint": "Side fraction >= unchanged immersed original at the same continuation stage"})
    np.savez_compressed(output / "mesh.npz", points=fem.points, triangles=fem.triangles,
                        velocity_points=fem.velocity_points, velocity_dofs=fem.velocity_dofs,
                        cell_index=domain.cell_index, density_x_mm=pipeline.grid.x_edges_mm,
                        density_y_mm=pipeline.grid.y_edges_mm)
    history, reports, figure_tasks = [], [], []
    figure_errors = []
    with (output / "progress.jsonl").open("x", encoding="utf-8") as journal, ThreadPoolExecutor(max_workers=1) as plots:
        for number, stage in enumerate(STAGES, start=1):
            folder = output / f"stage_{number:02d}_beta_{stage.beta:g}"
            design, ref_state, baseline, candidate, report = run_stage(
                folder, domain, pipeline, stage, number, design, ref_state, journal, history)
            reports.append(report)
            figures = folder / "images"
            figures.mkdir()
            figure_tasks.append(plots.submit(plot_stage, figures, domain, pipeline, baseline, candidate, tuple(history)))
        for number, task in enumerate(figure_tasks, start=1):
            try:
                task.result()
            except Exception as error:
                figure_errors.append({"stage": number, "error": f"{type(error).__name__}: {error}"})
    return {"sol_job_id": job, "stages": reports, "elapsed_seconds": perf_counter() - started,
            "figure_errors": figure_errors, "mma_unique_evaluations_total": sum(r["mma_unique_evaluations"] for r in reports),
            "mma_objective_calls_total": sum(r["mma_objective_calls"] for r in reports),
            "scope": "Exploratory immersed performance, with reconstructed geometry for inspection"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--local", action="store_true", help="Run the exploratory pilot on the local machine")
    args = parser.parse_args()
    if not args.local and not os.environ.get("SLURM_JOB_ID"):
        parser.error("Use --local or run under Slurm")
    args.output.mkdir(parents=True, exist_ok=False)
    completed = write_report(args.output, "summary", lambda: run(args.output))
    summary = json.loads((args.output / "summary.json").read_text(encoding="utf-8"))
    raise SystemExit(0 if completed and not summary["results"]["figure_errors"] else 1)


if __name__ == "__main__":
    main()
