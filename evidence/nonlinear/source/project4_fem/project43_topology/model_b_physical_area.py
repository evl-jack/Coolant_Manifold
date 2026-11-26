"""Run Model B with a conservative reconstruction-area constraint on ASU Sol."""

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
from .kkt_audit import audit_kkt
from .optimize import optimize
from .physical_area import ReconstructedAreaUpperBound
from .relaxed_flow import (
    BANDS, MMA_BUDGET, RESTORE_BUDGET, BandedEvaluator, candidate_designs, independent_replay, save, setup,
)


AREA_TRANSITION = .02
AREA_REGULARIZATION = 1e-6
AREA_ORDER = ("side_fraction_min", "side_fraction_max", "projected_area_mm2", "reconstruction_area_upper_mm2",
              "inlet_to_main", "inlet_to_side")


class AreaConstrainedEvaluator:
    """Keep the original hydraulic rows and the two hard restoration graph rows."""

    def __init__(self, hydraulics, graph, area):
        self.hydraulics = hydraulics
        self.connected = ConnectivityDiagnosticEvaluator(hydraulics, graph)
        self.area = area

    def evaluate(self, design):
        response = self.connected.evaluate(design)
        density = self.hydraulics.density(design)
        upper, projected_gradient = self.area.evaluate(density.projected)
        gradient = self.hydraulics.pipeline.filter_pullback(density.projection_derivative * projected_gradient)
        constraints = np.r_[response.constraints[:3], upper - MAX_AREA_MM2, response.constraints[3:]]
        jacobian = np.vstack((response.jacobian[:3], gradient, response.jacobian[3:]))
        return replace(response, constraints=constraints, jacobian=jacobian)


def common_reference(common, domain):
    immersed_path = common / "baseline.npz"
    fitted_path = common / "sharp_fitted_baseline.json"
    with np.load(immersed_path) as data:
        reference = float(data["power"])
        values = data["state"].copy()
    if not np.isfinite(reference) or reference <= 0. or values.shape != (domain.fem.size,) or not np.isfinite(values).all():
        raise ValueError("The common immersed reference does not match the original FEM")
    archive = json.loads(fitted_path.read_text(encoding="utf-8"))
    if not archive["completed"] or not archive["results"]["fitted_diagnostic_passed"]:
        raise ValueError("The common sharp fitted reference failed its mesh or state check")
    fitted_power = archive["results"]["levels"][-1]["mechanical_power_W_per_m"]
    if not np.isfinite(fitted_power) or fitted_power <= 0.:
        raise ValueError("The common sharp fitted reference power must be positive")
    free = np.setdiff1d(np.arange(domain.fem.size), domain.fixed)
    state = State(values, [], None, free)
    digests = {path.name: sha256(path.read_bytes()).hexdigest() for path in (immersed_path, fitted_path)}
    return reference, fitted_power, state, digests


def replay_candidate(optimization, original, graph, design, candidate_dir, reference_fitted, area):
    """Replay the standard physical gate and audit all six constraints using the original FEM."""
    report = independent_replay(original, graph, design, candidate_dir, reference_fitted)
    density = original.density(design)
    replay = AreaConstrainedEvaluator(original, graph, area)
    comparison = replay.evaluate(design)
    upper_area = comparison.constraints[3] + MAX_AREA_MM2
    exact_area = report["geometry"]["physical_fluid_area_mm2"]
    np.testing.assert_allclose(np.delete(comparison.constraints, 3), report["inequalities"], rtol=1e-10, atol=1e-9)
    roundoff = 64. * np.finfo(float).eps * max(1., abs(upper_area), abs(exact_area))
    if exact_area > upper_area + roundoff:
        raise RuntimeError(f"Exact reconstructed area {exact_area:.12g} exceeds the conservative bound {upper_area:.12g}")
    audit = audit_kkt(design, comparison.gradient, comparison.constraints, comparison.jacobian)
    accelerated = optimization.evaluate(design)
    report.update(
        projected_density=density.projected, filtered_density=density.filtered,
        conservative_area_upper_mm2=upper_area, conservative_area_margin_mm2=upper_area - exact_area,
        augmented_inequalities=comparison.constraints, augmented_kkt=asdict(audit),
        accelerated_replay={
            "objective": accelerated.objective, "inequalities": accelerated.constraints,
            "objective_difference": accelerated.objective - comparison.objective,
            "constraint_difference": accelerated.constraints - comparison.constraints,
        },
        handoff_eligible=(
            report["immersed_feasible"] and report["geometry"]["connected_ports"] and exact_area <= MAX_AREA_MM2
            and max(comparison.constraints) <= FEASIBILITY_TOLERANCE
        ),
    )
    return report


def run(output, common, max_stages):
    from Project43_Accelerated.hybrid import accelerated_domain

    domain, pipeline, graph = setup()
    if pipeline.grid.size != 240:
        raise RuntimeError("Model B requires the original 240-variable full-density grid")
    accelerated = accelerated_domain(domain, device="cpu")
    domain.fem.inertia = accelerated.fem.inertia = REYNOLDS
    reference, reference_fitted, seed, digests = common_reference(common, domain)
    incoming = pipeline.grid.base_raw.copy()
    area = ReconstructedAreaUpperBound(pipeline.grid, transition=AREA_TRANSITION, regularization=AREA_REGULARIZATION)
    planned = [Stage(.01, beta, alpha) for beta, alpha in ROUTES["coupled"]]
    stages = planned[:max_stages]
    band = BANDS["model-b"]
    save(output, "setup", {
        "model": "model-b", "diagnostic_only": True, "flow_band": band, "inertia": REYNOLDS,
        "design_cells": pipeline.grid.size, "design_spacing_mm": pipeline.grid.spacing_mm,
        "fem_spacing_mm": domain.max_spacing_mm, "starting_raw_sharp_t_design": incoming,
        "common_baseline": str(common.resolve()), "common_archive_sha256": digests,
        "normalization_W_per_m": reference, "fitted_sharp_power_W_per_m": reference_fitted,
        "area_cap_mm2": MAX_AREA_MM2, "area_transition": AREA_TRANSITION, "area_regularization": AREA_REGULARIZATION,
        "constraint_order": AREA_ORDER, "conservative_area_constraint_units": "mm²", "projection_threshold": .5,
        "mma_budget_per_stage": MMA_BUDGET, "restoration_budget_per_attempt": RESTORE_BUDGET,
        "planned_stages": [asdict(stage) for stage in planned], "stages": [asdict(stage) for stage in stages],
    })

    stage_reports = []
    stop_reason = "scheduled stages completed"
    for index, stage in enumerate(stages):
        folder = output / f"stage_{index:02d}"
        folder.mkdir()
        flow = BandedEvaluator(accelerated, pipeline, stage, reference, band=band, initial_state=seed)
        evaluator = AreaConstrainedEvaluator(flow, graph, area)
        start, initial, restoration = restore_stage(evaluator, flow, incoming, RESTORE_BUDGET, True)
        save(folder, "restoration", {"stage": asdict(stage), "initial": initial, "attempts": restoration,
                                     "accepted_design": start})
        if start is None:
            stop_reason = "no augmented-feasible start after restoration"
            break
        outcome, work = measured(lambda: optimize(evaluator, start, algorithm="mma", max_evaluations=MMA_BUDGET), flow)
        save(folder, "mma", {"stage": asdict(stage), "result": asdict(outcome), "cost": work})
        native_converged = outcome.optimizer_success and outcome.status in (1, 2, 3, 4)
        candidates = []
        for label, design in candidate_designs(outcome, start):
            candidate = folder / label
            candidate.mkdir()
            relative = candidate.relative_to(output).as_posix()
            original = BandedEvaluator(domain, pipeline, stage, reference, band=band, initial_state=seed)
            # A fresh accelerated replay leaves the audited optimization work counters intact.
            paired_flow = BandedEvaluator(accelerated, pipeline, stage, reference, band=band, initial_state=seed)
            paired = AreaConstrainedEvaluator(paired_flow, graph, area)
            try:
                report, cost = measured(
                    lambda: replay_candidate(paired, original, graph, design, candidate, reference_fitted, area), original)
            except (RuntimeError, ValueError, AssertionError, np.linalg.LinAlgError) as error:
                save(candidate, "checkpoint", {"stage": asdict(stage), "design": design,
                                               "error": f"{type(error).__name__}: {error}"})
                continue
            converged_physical = (label == "terminal" and native_converged and report["handoff_eligible"]
                                  and report["fitted_passed"] and report["augmented_kkt"]["verified"])
            report.update(stage=asdict(stage), audit_cost=cost, native_terminal_converged=native_converged,
                          converged_physical_candidate=converged_physical)
            save(candidate, "checkpoint", report)
            candidates.append((relative, report, original.state))

        handoffs = [item for item in candidates if item[1]["handoff_eligible"]]
        selected = min(handoffs, key=lambda item: (not item[1]["converged_physical_candidate"],
                       not item[1]["fitted_passed"], item[1]["objective"])) if handoffs else None
        summary = {"stage": asdict(stage), "mma_status": outcome.status, "mma_message": outcome.message,
                   "native_terminal_converged": native_converged, "continuation_handoff": None,
                   "fitted_screen_passes": [item[0] for item in candidates if item[1]["fitted_passed"]],
                   "converged_physical_candidate": None}
        if selected is not None:
            path, report, seed = selected
            incoming = report["design"]
            summary["continuation_handoff"] = path
            if report["converged_physical_candidate"]:
                summary["converged_physical_candidate"] = {
                    "checkpoint": path, "fitted_power_W_per_m": report["fitted_power_W_per_m"],
                    "fitted_power_ratio_vs_sharp": report["fitted_power_ratio_vs_sharp"],
                    "fitted_side_share": report["fitted"]["levels"][-1]["side_share"],
                    "fitted_area_mm2": report["geometry"]["physical_fluid_area_mm2"],
                }
        save(folder, "selected", summary)
        stage_reports.append(summary)
        if selected is None:
            stop_reason = "no independently replayed handoff passed the augmented and exact-geometry checks"
            break

    reached_final_stage = len(stage_reports) == len(planned)
    accepted = stage_reports[-1]["converged_physical_candidate"] if reached_final_stage else None
    return {"model": "model-b", "diagnostic_only": True, "stop_reason": stop_reason,
            "stage_results": stage_reports, "reached_final_stage": reached_final_stage, "accepted_endpoint": accepted}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="New per-job evidence directory")
    parser.add_argument("--common", type=Path, required=True, help="Existing shared baseline directory")
    parser.add_argument("--max-stages", type=int, help="Cap an exploratory run below the complete schedule")
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical diagnostic as a Sol Slurm job")
    if args.max_stages is not None and not 1 <= args.max_stages <= len(ROUTES["coupled"]):
        parser.error("Stage cap is outside the established continuation schedule")
    args.output.mkdir(parents=True, exist_ok=False)
    write_report(args.output, "environment", environment)
    if not write_report(args.output, "summary", lambda: run(args.output, args.common, args.max_stages)):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
