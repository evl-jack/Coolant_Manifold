"""Verify the reviewed Phase I and Phase II topology workflow on ASU Sol."""

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
from resource import RUSAGE_SELF, getrusage
from time import perf_counter, process_time

import numpy as np

from .density import DensityPipeline, DesignGrid, ReducedDensity
from .domain import build_domain
from .evaluator import Evaluator, Stage, baseline_reference
from .gate2 import _optimizer_result, environment, write_report
from .kkt_fixtures import run_kkt_fixtures
from .mma_diagnostics import CONSTRAINTS, EXPECTED_DESIGN, KnownProblem, independent_values
from .optimize import optimize
from .restoration import restore


PILOT_EVALUATIONS = 250
FEASIBILITY_TOLERANCE = 1e-6
COUNTERS = ("forward_attempts", "forward_solves", "adjoint_rhs", "newton_iterations")


def measured(function, evaluator=None):
    before = {item: getattr(evaluator, item, 0) for item in COUNTERS}
    wall_start, cpu_start = perf_counter(), process_time()
    result = function()
    cost = {item: getattr(evaluator, item, 0) - before[item] for item in COUNTERS}
    cost.update(wall_seconds=perf_counter() - wall_start, cpu_hours=(process_time() - cpu_start) / 3600,
                process_peak_rss_kib_at_end=getrusage(RUSAGE_SELF).ru_maxrss)
    return result, cost


def restore_start(evaluator, initial):
    reports = []
    for margin in (0.05, 0.0):
        outcome, cost = measured(lambda: restore(evaluator, initial, margin=margin,
                                                 max_evaluations=PILOT_EVALUATIONS), evaluator)
        reports.append(dict(result=asdict(outcome), cost=cost))
        if outcome.accepted_design is not None:
            return np.asarray(outcome.accepted_design), reports
    return None, reports


def phase_two(evaluator, start, algorithm):
    outcome, optimization_cost = measured(lambda: optimize(evaluator, start, algorithm=algorithm,
                                                            max_evaluations=PILOT_EVALUATIONS), evaluator)
    report, audit_cost = measured(lambda: _optimizer_result(evaluator, outcome), evaluator)
    report.update(optimization_cost=optimization_cost, audit_cost=audit_cost,
                  optimization_history=[asdict(item) for item in outcome.evaluation_history],
                  major_iteration_history=outcome.iteration_history)
    return report


def analytical():
    names = CONSTRAINTS["gate2"]
    accepted, restoration = restore_start(KnownProblem(names), np.zeros(2))
    phases = {}
    if accepted is not None:
        for algorithm in ("slsqp", "mma"):
            report = phase_two(KnownProblem(names), accepted, algorithm)
            design = np.asarray(report["terminal_design"])
            actual = independent_values(design, names)
            report.update(independent_objective=float(actual[0]), independent_inequalities=actual[1:].tolist(),
                          design_error=float(np.max(np.abs(design - EXPECTED_DESIGN))))
            report["passed"] = bool(report["design_error"] <= 1e-4 and report["kkt"] is not None
                                    and report["kkt"]["verified"] and max(actual[1:]) <= FEASIBILITY_TOLERANCE)
            phases[algorithm] = report
    return dict(initial=[0., 0.], restoration=restoration, phase_two=phases,
                passed=bool(phases and all(item["passed"] for item in phases.values())))


def topology(output, mode):
    grid = DesignGrid(0.5)
    pipeline = DensityPipeline(grid)
    domain = build_domain(grid, 0.5)
    stage = Stage(0.01, 0., 1e5)
    (power_reference, seed), setup_cost = measured(lambda: baseline_reference(domain, pipeline, stage.alpha_max))
    setup_cost.update(forward_solves=3, forward_attempts=3, newton_iterations=None)
    reduced = ReducedDensity(grid) if mode == "reduced" else None
    initial = np.zeros(20) if reduced is not None else grid.base_raw
    evaluator = Evaluator(domain, pipeline, stage, power_reference, reduced=reduced, initial_state=seed)
    accepted, restoration = restore_start(evaluator, initial)
    result = dict(mode=mode, initial=initial.tolist(), restoration=restoration, setup_cost=setup_cost,
                  mesh=dict(triangles=len(domain.fem.triangles), dofs=domain.fem.size, design_cells=grid.size),
                  stage=asdict(stage), power_reference_w_per_m=power_reference)
    write_report(output, f"{mode}_phase_one", lambda: result)
    if accepted is None:
        return dict(restored=False)

    audit = Evaluator(domain, pipeline, stage, power_reference, reduced=reduced, initial_state=seed)
    anchor, audit_cost = measured(lambda: audit.evaluate(accepted), audit)
    summary = dict(restored=bool(max(anchor.constraints) <= FEASIBILITY_TOLERANCE), anchor=accepted.tolist(),
                   anchor_objective=anchor.objective, anchor_constraints=anchor.constraints.tolist(),
                   anchor_audit_cost=audit_cost)
    if not summary["restored"]:
        return summary

    if reduced is not None:
        equivalent = Evaluator(domain, pipeline, stage, power_reference, initial_state=audit.state)
        equivalent_result, transfer_cost = measured(lambda: equivalent.evaluate(reduced.raw(accepted)), equivalent)
        summary.update(full_space_transfer=dict(raw=reduced.raw(accepted).tolist(), cost=transfer_cost,
                            state_max_difference=float(np.max(np.abs(equivalent.state.values - audit.state.values))),
                            objective_difference=abs(equivalent_result.objective-anchor.objective),
                            constraint_max_difference=float(np.max(np.abs(equivalent_result.constraints-anchor.constraints)))))
    else:
        for algorithm in ("slsqp", "mma"):
            branch = Evaluator(domain, pipeline, stage, power_reference, initial_state=audit.state)
            write_report(output, f"phase_two_{algorithm}", lambda: phase_two(branch, accepted, algorithm))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("mode", choices=("analytical", "full", "reduced"))
    arguments = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical verification as a Sol Slurm job")
    arguments.output.mkdir(parents=True)
    write_report(arguments.output, "environment", environment)
    if arguments.mode == "analytical":
        write_report(arguments.output, "analytical", analytical)
        write_report(arguments.output, "kkt", run_kkt_fixtures)
        for name in ("analytical", "kkt"):
            report = json.loads((arguments.output / f"{name}.json").read_text(encoding="utf-8"))
            if not report["completed"] or not report["results"]["passed"]:
                raise SystemExit(1)
    else:
        write_report(arguments.output, "restoration_summary", lambda: topology(arguments.output, arguments.mode))


if __name__ == "__main__":
    main()
