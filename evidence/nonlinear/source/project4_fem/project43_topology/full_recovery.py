"""Replay Phase II best feasible checkpoints and recover 240-cell endpoints on Sol."""

import argparse
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path

import numpy as np

from .density import DensityPipeline, DesignGrid
from .domain import build_domain
from .evaluator import Evaluator, Stage, baseline_reference
from .gate2 import environment, write_report
from .gate2_restoration import FEASIBILITY_TOLERANCE, measured
from .kkt_audit import audit_kkt
from .optimize import optimize
from .verify_flow import STEPS, _direction_check, _power_and_flux


PER_RESTART_EVALUATIONS = 250
MAX_RESTARTS = 3
EXPERIMENT_STAGE = Stage(0.01, 0.0, 1e5)
STATIONARITY_BRACKET_WIDTH = 1e-8
STATIONARITY_BISECTIONS = 20


def read_archive(directory, name):
    path = directory / f"{name}.json"
    payload = path.read_bytes()
    data = json.loads(payload)
    if not data["completed"]:
        raise ValueError(f"Archived experiment did not complete: {path}")
    return data["results"], {"file": str(path), "sha256": hashlib.sha256(payload).hexdigest()}


def archived_checkpoint(archive, size):
    design = np.asarray(archive["best_feasible_design"], dtype=float)
    if design.shape != (size,) or not np.all(np.isfinite(design)) or np.any((design < 0) | (design > 1)):
        raise ValueError("Archive does not contain a finite in-box best feasible checkpoint")
    feasible = [row for row in archive["optimization_history"]
                if max(row["inequalities"]) <= FEASIBILITY_TOLERANCE
                and all(0 <= value <= 1 for value in row["design"])]
    best = min(feasible, key=lambda row: row["objective"], default=None)
    if best is None or not np.array_equal(design, np.asarray(best["design"])):
        raise ValueError("Archived checkpoint does not match the lowest-objective feasible history row")
    if not 1 <= archive["unique_evaluations"] <= PER_RESTART_EVALUATIONS:
        raise ValueError("Original Phase II evaluation count is outside its approved budget")
    return design, best


def stationarity_bracket(design, gradient, inequalities, jacobian):
    """Numerically bracket the best stationarity with the original multiplier conditions."""
    lower = 0.0
    # Zero multipliers satisfy dual feasibility and complementarity exactly.
    upper = float(np.max(np.abs(design - np.clip(design - gradient, 0.0, 1.0))))
    upper_multipliers = np.zeros(inequalities.size)
    probes = []
    stop_reason = "probe limit"
    for _ in range(STATIONARITY_BISECTIONS):
        if upper - lower <= STATIONARITY_BRACKET_WIDTH:
            stop_reason = "target bracket width"
            break
        midpoint = 0.5 * (lower + upper)
        audit = audit_kkt(design, gradient, inequalities, jacobian, stationarity_tolerance=midpoint)
        probes.append({"threshold": midpoint, "lp_status": audit.lp_status,
                       "witness_stationarity": audit.stationarity_inf, "verified": audit.verified})
        if audit.verified:
            upper = audit.stationarity_inf
            upper_multipliers = audit.multipliers
            if upper < lower:
                lower = 0.0
                stop_reason = "witness contradicted numerical LP lower bound"
                break
        elif audit.lp_status == 2:
            lower = midpoint
        else:
            # A numerical LP status without a passing direct residual does not
            # establish feasibility or infeasibility of its exact inequalities.
            stop_reason = "indeterminate numerical LP or residual"
            break
    return {"lower_from_lp_infeasibility": lower, "upper_from_verified_witness": upper,
            "upper_multipliers": upper_multipliers, "bracket_width": upper-lower,
            "probes": probes, "stop_reason": stop_reason, "informational_only": True}


def audit_checkpoint(evaluator, design):
    value = evaluator.evaluate(design)
    power, flux = _power_and_flux(evaluator.domain, evaluator.state.values)
    density = evaluator.density(design)
    area = 41.0 + evaluator.pipeline.grid.cell_area_mm2 * float(np.sum(density.projected))
    split = float(flux[2] / flux[0])
    objective = power / evaluator.power_reference
    inequalities = np.array([(0.499-split)/0.001, (split-0.501)/0.001, area-83.0])
    kkt = audit_kkt(design, value.gradient, inequalities, value.jacobian)
    return {"design": design, "objective": objective, "inequalities": inequalities,
            "power_w_per_m": power, "side_fraction": split, "area_mm2": area,
            "gray_measure": value.gray_measure, "mass_imbalance": value.mass_imbalance,
            "nonlinear_residual": value.nonlinear_residual,
            "engineering_feasible": kkt.primal_violation <= FEASIBILITY_TOLERANCE,
            "response_replay_max_difference": float(max(abs(objective-value.objective),
                np.max(np.abs(inequalities-value.constraints)))),
            "kkt": asdict(kkt),
            "stationarity_bracket": (stationarity_bracket(design, value.gradient, inequalities, value.jacobian)
                                     if kkt.primal_violation <= FEASIBILITY_TOLERANCE else None)}


def derivative_probe(evaluator, design):
    rng = np.random.default_rng(432395874)
    random = rng.uniform(0.1, 1.0, size=design.size)
    inward = random * np.where(design <= 0.5, 1.0, -1.0)
    directions = {"inward": (inward, True)}
    eligible = (design >= STEPS[0]) & (design <= 1.0-STEPS[0])
    if np.any(eligible):
        centered = np.where(eligible, random * rng.choice([-1.0, 1.0], size=design.size), 0.0)
        directions["interior_supported"] = (centered, False)
    return {name: {"direction": direction, **_direction_check(evaluator, design, direction, at_bound=at_bound)}
            for name, (direction, at_bound) in directions.items()}


def replay(evaluator_factory, design, *, probe):
    evaluator = evaluator_factory()
    audit, audit_cost = measured(lambda: audit_checkpoint(evaluator, design), evaluator)
    audit["audit_cost"] = audit_cost
    if probe and audit["engineering_feasible"]:
        attach_probe(audit, evaluator)
    return audit, evaluator


def attach_probe(audit, evaluator):
    def check():
        try:
            return derivative_probe(evaluator, audit["design"])
        except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
            return {"error": f"{type(error).__name__}: {error}"}

    probes, verification_cost = measured(check, evaluator)
    audit["derivative_probes"] = probes
    audit["derivative_verification_cost"] = verification_cost


def recover_branch(algorithm, archive_dir, output, evaluator_factory, size):
    archive, archive_file = read_archive(archive_dir, f"phase_two_{algorithm}")
    start, best = archived_checkpoint(archive, size)
    current, evaluator = replay(evaluator_factory, start, probe=True)
    checkpoint_file = f"{algorithm}_archive_checkpoint"
    if not write_report(output, checkpoint_file, lambda: {
        "archive": archive_file, "archived_best_trial": best, "archived_terminal_design": archive["terminal_design"],
        "archived_terminal_kkt": archive["kkt"], "archived_optimization_cost": archive["optimization_cost"],
        "checkpoint": current,
    }):
        raise RuntimeError("Could not save the archive checkpoint audit")

    evaluations = archive["unique_evaluations"]
    reports = [f"{checkpoint_file}.json"]
    reason = "archived best trial failed independent feasibility"
    unchanged_restarts = 0
    for restart in range(1, MAX_RESTARTS + 1):
        if not current["engineering_feasible"] or current["kkt"]["verified"]:
            break

        # Recreate the optimizer and seed the FEM state from the audited design.
        iteration_evaluator = evaluator_factory(evaluator.state)
        outcome, cost = measured(lambda: optimize(iteration_evaluator, current["design"], algorithm=algorithm,
                                                  max_evaluations=PER_RESTART_EVALUATIONS), iteration_evaluator)
        evaluations += outcome.unique_evaluations
        views = {}
        for name, design in (("terminal", outcome.terminal_design), ("best_feasible_trial", outcome.best_feasible_design)):
            if design is None or not np.all(np.isfinite(design)) or np.any((design < 0) | (design > 1)):
                continue
            if np.array_equal(design, current["design"]):
                views[name] = current, evaluator
            elif name == "best_feasible_trial" and "terminal" in views and np.array_equal(design, outcome.terminal_design):
                views[name] = views["terminal"]
            else:
                try:
                    views[name] = replay(evaluator_factory, design, probe=False)
                except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
                    views[name] = {"design": design, "engineering_feasible": False,
                                   "error": f"{type(error).__name__}: {error}"}, None

        candidates = {"incoming_checkpoint": (current, evaluator),
                      **{name: item for name, item in views.items() if item[0]["engineering_feasible"]}}
        stationary = [name for name, (report, _) in candidates.items() if report["kkt"]["verified"]]
        selected = min(stationary or candidates, key=lambda name: candidates[name][0]["objective"])
        candidate, next_evaluator = candidates[selected]
        if "derivative_probes" not in candidate:
            attach_probe(candidate, next_evaluator)
        if outcome.failed_design is None and np.array_equal(candidate["design"], current["design"]):
            unchanged_restarts += 1
        else:
            unchanged_restarts = 0
        restart_record = {
            "algorithm": algorithm, "restart": restart, "evaluation_budget": PER_RESTART_EVALUATIONS,
            "cumulative_optimizer_evaluations": evaluations,
            "native_state_preserved": False, "warm_start_fem_from_incoming_checkpoint": True,
            "incoming_design": current["design"], "incoming_objective": current["objective"],
            "optimization_cost": cost, "optimizer": asdict(outcome),
            "independent_audits": {name: report for name, (report, _) in views.items()},
            "selected_checkpoint_source": selected, "selected_checkpoint": candidate,
            "consecutive_restarts_without_checkpoint_change": unchanged_restarts,
        }
        report_name = f"{algorithm}_restart_{restart:02}"
        if not write_report(output, report_name, lambda: restart_record):
            raise RuntimeError(f"Could not save restart {restart} evidence")
        reports.append(f"{report_name}.json")
        current, evaluator = candidate, next_evaluator
        reason = "restart evaluation budget exhausted"
        if unchanged_restarts >= 2:
            reason = "two fresh restarts left the exact same independently audited checkpoint"
            break

    if current["kkt"]["verified"]:
        reason = "original KKT thresholds verified"
    return {"algorithm": algorithm, "archive": archive_file, "reports": reports,
            "optimizer_evaluations_including_archive": evaluations,
            "native_state_preserved": False, "stop_reason": reason, "final_checkpoint": current,
            "passed": current["kkt"]["verified"]}


def run(archive_dir, output, algorithms, endpoint_root=None):
    grid = DesignGrid(0.5)
    pipeline = DensityPipeline(grid)
    domain = build_domain(grid, 0.5)
    archived_setup, setup_file = read_archive(archive_dir, "full_phase_one")
    mesh = {"triangles": len(domain.fem.triangles), "dofs": domain.fem.size, "design_cells": grid.size}
    if mesh != archived_setup["mesh"] or asdict(EXPERIMENT_STAGE) != archived_setup["stage"]:
        raise ValueError("Recovery mesh and stage must exactly match the archived full-topology fixture")

    (reference, seed), setup_cost = measured(lambda: baseline_reference(domain, pipeline, EXPERIMENT_STAGE.alpha_max))
    setup_cost.update(forward_solves=3, forward_attempts=3, newton_iterations=None)
    archived_reference = archived_setup["power_reference_w_per_m"]
    if not np.isclose(reference, archived_reference, rtol=1e-10, atol=0.0):
        raise ValueError("Independent baseline power does not match the fixed archived P0")
    if not write_report(output, "setup", lambda: {"archive": setup_file, "mesh": mesh,
        "stage": asdict(EXPERIMENT_STAGE), "archived_power_reference_w_per_m": archived_reference,
        "replayed_power_reference_w_per_m": reference, "setup_cost": setup_cost}):
        raise RuntimeError("Could not save matched setup evidence")

    def evaluator_factory(initial_state=seed):
        return Evaluator(domain, pipeline, EXPERIMENT_STAGE, archived_reference, initial_state=initial_state)

    results = {}
    for algorithm in algorithms:
        name = f"{algorithm}_summary"
        if endpoint_root is None:
            write_report(output, name,
                         lambda: recover_branch(algorithm, archive_dir, output, evaluator_factory, grid.size))
        else:
            def endpoint_replay():
                report, archive_path = read_archive(endpoint_root / algorithm, name)
                design = np.asarray(report["final_checkpoint"]["design"], dtype=float)
                audit, _ = replay(evaluator_factory, design, probe=True)
                probes = audit.get("derivative_probes", {})
                return {"archive": archive_path, "final_checkpoint": audit,
                        "passed": bool(probes and all(isinstance(value, dict) and value.get("passed", False)
                                                     for value in probes.values()))}

            write_report(output, name, endpoint_replay)
        results[algorithm] = json.loads((output / f"{name}.json").read_text(encoding="utf-8"))
    return {"branches": results, "job_id": os.environ["SLURM_JOB_ID"], "archive": str(archive_dir),
            "all_branches_completed": all(result["completed"] for result in results.values()),
            "all_branches_verified": all(result["completed"] and result["results"]["passed"]
                                         for result in results.values())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="Archived 2395874 full-topology report directory")
    parser.add_argument("output", type=Path, help="New recovery evidence directory")
    parser.add_argument("--algorithm", choices=("slsqp", "mma"))
    parser.add_argument("--endpoint-root", type=Path, help="Replay derivatives at archived final checkpoints only")
    arguments = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical verification as a Sol Slurm job")
    arguments.output.mkdir(parents=True)
    write_report(arguments.output, "environment", environment)
    algorithms = (arguments.algorithm,) if arguments.algorithm else ("slsqp", "mma")
    completed = write_report(arguments.output, "recovery_summary",
                             lambda: run(arguments.archive, arguments.output, algorithms, arguments.endpoint_root))
    summary = json.loads((arguments.output / "recovery_summary.json").read_text(encoding="utf-8"))
    if not completed or not summary["results"]["all_branches_completed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
