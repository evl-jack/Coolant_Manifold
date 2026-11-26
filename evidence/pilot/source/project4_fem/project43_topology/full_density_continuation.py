"""Bounded full-density projection and Brinkman diagnostics on ASU Sol."""

import argparse
from dataclasses import asdict
import os
from pathlib import Path

import numpy as np

from .density import DesignGrid, DensityPipeline
from .domain import build_domain
from .evaluator import Evaluator, Stage, baseline_reference
from .full_recovery import audit_checkpoint
from .gate2 import environment, write_report
from .gate2_restoration import FEASIBILITY_TOLERANCE, measured
from .kkt_audit import audit_kkt
from .optimize import optimize
from .reconstruction import reconstruct, solve_fitted, write_reconstruction
from .restoration import restore


ROUTES = {
    "projection": ((0., 1e5), (1., 1e5), (2., 1e5), (4., 1e5), (8., 1e5)),
    "brinkman": ((0., 1e5), (0., 3e5), (0., 1e6)),
    "coupled": ((0., 1e5), (1., 1e5), (1., 3e5), (2., 3e5), (2., 1e6), (4., 1e6), (8., 1e6)),
}
FITTED_SPACINGS_MM = (0.5, 0.25, 0.125)
FLOW_BAND = (0.499, 0.501)
MAX_AREA_MM2 = 83.
FITTED_SPLIT_CHANGE_PP = 0.02
FITTED_POWER_RELATIVE_CHANGE = 0.005


def fitted_screen(geometry, flow_band=FLOW_BAND, spacings=FITTED_SPACINGS_MM, *, area_cap_mm2=MAX_AREA_MM2):
    """Independently check the hard-wall geometry on three fixed mesh levels."""
    if geometry.port_component is None:
        return {"attempted": False, "reason": "inlet and both outlets are not connected"}
    if area_cap_mm2 is not None and geometry.area_mm2 > area_cap_mm2:
        return {"attempted": False, "reason": "reconstructed fluid area exceeds the cap"}

    replays = []
    for spacing in spacings:
        try:
            (response, _, _), cost = measured(lambda: solve_fitted(geometry, spacing))
        except Exception as error:
            return {"attempted": True, "completed": False, "spacing_mm": spacing,
                    "error": f"{type(error).__name__}: {error}", "levels": replays}
        replays.append(dict(response, in_physical_flow_band=flow_band[0] <= response["side_share"] <= flow_band[1],
                            cost=cost))
    return summarize_fitted_replays(replays, flow_band)


def summarize_fitted_replays(replays, flow_band):
    """Apply the existing mesh and flow acceptance tests to the ordered levels."""
    changes = []
    for coarse, fine in zip(replays[:-1], replays[1:]):
        split_pp = 100 * abs(fine["side_share"] - coarse["side_share"])
        coarse_power = abs(coarse["mechanical_power_W_per_m"])
        fine_power = abs(fine["mechanical_power_W_per_m"])
        if min(coarse_power, fine_power) <= 0.:
            power_change = None
        else:
            power_change = abs(fine["mechanical_power_W_per_m"] - coarse["mechanical_power_W_per_m"]) / fine_power
        changes.append({"flow_change_percentage_points": split_pp, "relative_power_change": power_change,
                        "passed": split_pp <= FITTED_SPLIT_CHANGE_PP and power_change is not None
                        and power_change <= FITTED_POWER_RELATIVE_CHANGE})
    split_uncertainty = max(change["flow_change_percentage_points"] for change in changes) / 100.
    finest_split = replays[-1]["side_share"]
    split_interval = (finest_split - split_uncertainty, finest_split + split_uncertainty)
    interval_in_band = flow_band[0] <= split_interval[0] and split_interval[1] <= flow_band[1]
    return {"attempted": True, "completed": True, "levels": replays, "successive_changes": changes,
            "conservative_split_uncertainty_fraction": split_uncertainty,
            "conservative_fitted_split_interval": split_interval,
            "conservative_split_interval_in_flow_band": interval_in_band,
            "fitted_diagnostic_passed": all(level["state_passed"] and level["in_physical_flow_band"] for level in replays)
            and all(change["passed"] for change in changes) and interval_in_band}


def optimizer_evaluator(hydraulics, connectivity_diagnostic):
    if connectivity_diagnostic:
        from .connectivity import ConnectivityDiagnosticEvaluator

        return ConnectivityDiagnosticEvaluator(hydraulics)
    return hydraulics


def evaluate_candidate(domain, pipeline, stage, reference, seed, design, output, evaluate_fitted,
                       connectivity_diagnostic):
    """Audit a fresh immersed replay and retain the independent contour evidence."""
    evaluator = Evaluator(domain, pipeline, stage, reference, initial_state=seed)
    checkpoint, cost = measured(lambda: audit_checkpoint(evaluator, design), evaluator)
    density = evaluator.density(design)
    checkpoint.update(audit_cost=cost, filtered_density=density.filtered, projected_density=density.projected,
                      resistance=density.resistance)
    if connectivity_diagnostic:
        augmented = optimizer_evaluator(evaluator, True)
        response = augmented.evaluate(design)
        certificate, _ = augmented.certificate(design)
        normalized = response.constraints[3:]
        inequalities = np.r_[checkpoint["inequalities"], normalized]
        graph_audit = audit_kkt(design, response.gradient, inequalities, response.jacobian)
        checkpoint["connectivity_diagnostic"] = {
            "port_order": ("main", "side"), "conductance": certificate.conductance,
            "normalized_inequalities": normalized, "strict_threshold_paths": certificate.hard_paths,
            "minimum_conductance": augmented.graph.minimum_conductance, "floor_budget": augmented.graph.floor_budget,
            "strict_threshold": augmented.graph.cutoff, "seam_counts_inlet_main_side": augmented.graph.seam_counts,
            "augmented_inequalities": inequalities, "augmented_kkt": asdict(graph_audit),
            "connectivity_only_satisfied": bool(max(normalized) <= FEASIBILITY_TOLERANCE
                                                and all(certificate.hard_paths)),
        }
    try:
        geometry = reconstruct(pipeline.grid, density.projected)
        checkpoint["reconstruction"] = write_reconstruction(geometry, output / "geometry.json")
    except Exception as error:
        checkpoint["reconstruction"] = {"error": f"{type(error).__name__}: {error}"}
    else:
        if evaluate_fitted:
            checkpoint["fitted_screen"] = fitted_screen(geometry)
    return checkpoint, evaluator.state


def restore_stage(evaluator, hydraulics, incoming, evaluation_budget, connectivity_diagnostic):
    """Change the warm start only when the new stage needs feasibility restoration."""
    value, initial_cost = measured(lambda: evaluator.evaluate(incoming), hydraulics)
    initial = {"design": incoming, "inequalities": value.constraints, "flow_fraction": value.side_fraction,
               "area_mm2": value.area_mm2, "cost": initial_cost}
    attempts = []
    if max(value.constraints) <= FEASIBILITY_TOLERANCE:
        return incoming, initial, attempts
    for margin in (0.05, 0.):
        outcome, cost = measured(lambda: restore(evaluator, incoming, margin=margin,
                                                 max_evaluations=evaluation_budget,
                                                 hard_tail=2 if connectivity_diagnostic else 0), hydraulics)
        attempts.append({"result": asdict(outcome), "cost": cost})
        if outcome.accepted_design is not None:
            return np.asarray(outcome.accepted_design), initial, attempts
    return None, initial, attempts


def run_stage(domain, pipeline, stage, reference, seed, incoming, output, mma_budget, restore_budget,
              connectivity_diagnostic):
    hydraulics = Evaluator(domain, pipeline, stage, reference, initial_state=seed)
    evaluator = optimizer_evaluator(hydraulics, connectivity_diagnostic)
    start, initial, attempts = restore_stage(evaluator, hydraulics, incoming, restore_budget,
                                            connectivity_diagnostic)
    if not write_report(output, "restoration", lambda: {"initial": initial, "attempts": attempts,
                                                        "accepted_design": start}):
        raise RuntimeError("Could not preserve stage restoration evidence")

    candidates = [("incoming", incoming)]
    if start is not None:
        candidates.append(("restored_start", start))
        outcome, cost = measured(lambda: optimize(evaluator, start, algorithm="mma", max_evaluations=mma_budget),
                                 hydraulics)
        if not write_report(output, "mma", lambda: {"result": asdict(outcome), "cost": cost}):
            raise RuntimeError("Could not preserve MMA history")
        candidates.extend((("terminal", outcome.terminal_design), ("best_feasible_trial", outcome.best_feasible_design)))

    audits = {}
    for name, design in reversed(candidates):
        if design is None or not np.isfinite(design).all() or np.any((design < 0.) | (design > 1.)):
            continue
        if any(np.array_equal(design, candidate) for candidate, _ in audits.values()):
            continue
        candidate_dir = output / name
        candidate_dir.mkdir()
        try:
            checkpoint, audit_seed = evaluate_candidate(
                domain, pipeline, stage, reference, seed, design, candidate_dir,
                evaluate_fitted=name in ("terminal", "best_feasible_trial"),
                connectivity_diagnostic=connectivity_diagnostic,
            )
        except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
            checkpoint, audit_seed = {"design": design, "engineering_feasible": False,
                                      "error": f"{type(error).__name__}: {error}"}, None
        audits[name] = (design, (checkpoint, audit_seed))
        if not write_report(candidate_dir, "checkpoint", lambda: checkpoint):
            raise RuntimeError(f"Could not preserve checkpoint {name}")

    def acceptable(report):
        return report["engineering_feasible"] and (not connectivity_diagnostic or (
            report.get("connectivity_diagnostic", {}).get("connectivity_only_satisfied", False)
            and report.get("reconstruction", {}).get("connected_ports", False)))

    def stationary_for_optimization(report):
        if connectivity_diagnostic:
            return report["connectivity_diagnostic"]["augmented_kkt"]["verified"]
        return report["kkt"]["verified"]

    feasible = {name: pair for name, pair in audits.items() if acceptable(pair[1][0])}
    stationary = {name: pair for name, pair in feasible.items() if stationary_for_optimization(pair[1][0])}
    eligible = stationary or feasible
    next_name = min(eligible, key=lambda name: eligible[name][1][0]["objective"]) if eligible else None
    snapshot = {"stage": asdict(stage), "selected_continuation_source": next_name,
                "candidates": {name: {"immersed_feasible": report["engineering_feasible"],
                                      "immersed_kkt_verified": report.get("kkt", {}).get("verified", False),
                                      "augmented_connectivity_kkt_verified": report.get("connectivity_diagnostic", {}).get(
                                          "augmented_kkt", {}).get("verified"),
                                      "strict_threshold_paths": report.get("connectivity_diagnostic", {}).get("strict_threshold_paths"),
                                      "side_fraction": report.get("side_fraction"),
                                      "area_mm2": report.get("area_mm2"),
                                      "connected_ports": report.get("reconstruction", {}).get("connected_ports"),
                                      "reconstructed_area_mm2": report.get("reconstruction", {}).get("physical_fluid_area_mm2"),
                                      "fitted_diagnostic_passed": report.get("fitted_screen", {}).get("fitted_diagnostic_passed")}
                               for name, (_, (report, _)) in audits.items()}}
    if not write_report(output, "summary", lambda: snapshot):
        raise RuntimeError("Could not preserve stage summary")
    if next_name is None:
        return snapshot, None, None
    next_design, (_, next_seed) = eligible[next_name]
    return snapshot, next_design, next_seed


def run(output, route, envelope, design_spacing, fem_spacing, mma_budget, restore_budget, connectivity_diagnostic):
    if envelope == "expanded":
        from .expanded_envelope import ExpandedDesignGrid, build_expanded_domain

        grid = ExpandedDesignGrid(design_spacing)
        domain = build_expanded_domain(grid, fem_spacing)
    else:
        grid = DesignGrid(design_spacing)
        domain = build_domain(grid, fem_spacing)
    pipeline = DensityPipeline(grid)
    reference, seed = baseline_reference(domain, pipeline, 1e5)
    incoming = grid.base_raw.copy()
    stages = [Stage(0.01, beta, alpha_max) for beta, alpha_max in ROUTES[route]]
    provenance = {"slurm_job_id": os.environ["SLURM_JOB_ID"], "route": route, "envelope": envelope,
                  "design_cells": grid.size, "design_spacing_mm": design_spacing, "fem_spacing_mm": fem_spacing,
                  "mesh_triangles": len(domain.fem.triangles), "fem_dofs": domain.fem.size,
                  "design_bounds_mm": [[grid.x_edges_mm[0], grid.x_edges_mm[-1]],
                                       [grid.y_edges_mm[0], grid.y_edges_mm[-1]]],
                  "fixed_reference_power_W_per_m": reference, "initial_raw_sharp_t_design": incoming,
                  "stages": [asdict(stage) for stage in stages], "mma_evaluations_per_stage": mma_budget,
                  "restoration_evaluations_per_attempt": restore_budget, "projection_threshold": pipeline.threshold,
                  "diagnostic_only": True, "connectivity_in_optimizer": connectivity_diagnostic,
                  "connectivity_excluded_from_restoration_slack": connectivity_diagnostic,
                  "displacement_correction_in_optimizer": False}
    if not write_report(output, "setup", lambda: provenance):
        raise RuntimeError("Could not preserve continuation setup")
    reports = []
    for number, stage in enumerate(stages):
        stage_dir = output / f"stage_{number:02d}"
        stage_dir.mkdir()
        snapshot, incoming, seed = run_stage(domain, pipeline, stage, reference, seed, incoming,
                                            stage_dir, mma_budget, restore_budget, connectivity_diagnostic)
        reports.append(snapshot)
        if incoming is None:
            break
    return {"reached_all_scheduled_stages": len(reports) == len(stages),
            "stop_reason": "no verified feasible handoff" if incoming is None else "diagnostic schedule complete",
            "diagnostic_stage_results": reports, "last_continuation_design": incoming}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="New per-job evidence directory")
    parser.add_argument("--route", choices=ROUTES, default="projection")
    parser.add_argument("--envelope", choices=("original", "expanded"), default="original")
    parser.add_argument("--design-spacing", type=float, default=0.5)
    parser.add_argument("--fem-spacing", type=float, default=0.5)
    parser.add_argument("--mma-evaluations", type=int, default=250)
    parser.add_argument("--restoration-evaluations", type=int, default=250)
    parser.add_argument("--connectivity-diagnostic", action="store_true",
                        help="Append experimental strict-threshold graph constraints")
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical diagnostic as a Sol Slurm job")
    if min(args.mma_evaluations, args.restoration_evaluations) < 1:
        parser.error("Evaluation budgets must be positive")
    args.output.mkdir(parents=True, exist_ok=False)
    write_report(args.output, "environment", environment)
    ok = write_report(args.output, "continuation_summary", lambda: run(
        args.output, args.route, args.envelope, args.design_spacing, args.fem_spacing,
        args.mma_evaluations, args.restoration_evaluations, args.connectivity_diagnostic))
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
