"""Run the staged, anchored wall-displacement optimization on ASU Sol."""

import argparse
from dataclasses import asdict
from hashlib import sha256
import json
import os
from pathlib import Path
from time import perf_counter, process_time

import numpy as np

from .density import DesignGrid
from .domain import build_domain
from .evaluator import REYNOLDS, Stage
from .gate2 import environment, write_report
from .gate2_restoration import COUNTERS, FEASIBILITY_TOLERANCE
from .kkt_audit import audit_kkt
from .optimize import optimize
from .reconstruction import write_reconstruction
from .relaxed_flow import save
from .restoration import restore
from .verify_flow import _power_and_flux
from .wall_displacement_density_v1 import WallDisplacementDensityV1
from .wall_displacement_fem_v1 import PDE_MODES, WallDisplacementEvaluatorV1
from .wall_displacement_geometry_v1 import audit_wall_displacement, reconstruct_clipped


STAGES = tuple(Stage(q, beta, 1e5) for q, beta in (
    (.01, 0.), (.05, 0.), (.10, 0.), (.10, 1.), (.10, 2.), (.10, 4.), (.10, 8.),
))
SIDE_BAND = (.10, .90)
SPLIT_SCALE = .001
RESTORE_BUDGET = 250
MMA_BUDGET = 500


class WorkLog:
    """Track successful and failed attempts across independent evaluator instances."""

    def __init__(self):
        self.entries = []

    def measure(self, label, evaluator, operation):
        before = {name: getattr(evaluator, name) for name in COUNTERS}
        wall, cpu = perf_counter(), process_time()
        try:
            return operation()
        finally:
            cost = {name: getattr(evaluator, name) - before[name] for name in COUNTERS}
            self.entries.append(dict(cost, operation=label, wall_seconds=perf_counter() - wall,
                                     cpu_hours=(process_time() - cpu) / 3600.))

    def totals(self):
        return {key: sum(entry[key] for entry in self.entries)
                for key in (*COUNTERS, "wall_seconds", "cpu_hours")}


def sharp_reference(common):
    path = common / "sharp_fitted_baseline.json"
    payload = path.read_bytes()
    archive = json.loads(payload)
    record = archive["results"] if archive["completed"] else {}
    if not record.get("fitted_diagnostic_passed"):
        raise ValueError("The sharp fitted-wall reference has not passed its existing verification checks")
    finest = record["levels"][-1]
    power = finest["mechanical_power_W_per_m"]
    if not np.isfinite(power) or power <= 0.:
        raise ValueError("The verified sharp fitted reference power must be finite and positive")
    return power, {"archive": str(path.resolve()), "sha256": sha256(payload).hexdigest(), "finest_level": finest}


def record_design(evaluator, design, output, *, stationarity, geometry_constraint=None):
    """Independently integrate the port quantities and retain the material fields."""
    value = evaluator.evaluate(design)
    density = evaluator.density(design)
    power, flux = _power_and_flux(evaluator.domain, evaluator.state.values)
    if not np.isfinite(flux).all() or flux[0] <= 0.:
        raise ValueError("The independent port integration did not give a positive finite inlet flux")
    fraction = float(flux[2] / flux[0])
    objective = power / evaluator.power_reference
    inequalities = np.array([(SIDE_BAND[0] - fraction) / SPLIT_SCALE, (fraction - SIDE_BAND[1]) / SPLIT_SCALE])
    if geometry_constraint is not None:
        inequalities = np.r_[inequalities, geometry_constraint.value(design)]
    np.testing.assert_allclose(value.flux, flux, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(value.objective, objective, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(value.constraints, inequalities, rtol=1e-10, atol=1e-9)
    np.savez_compressed(output / "fields.npz", compact_design=design, raw=evaluator.pipeline.raw(design),
                        filtered=density.filtered, projected=density.projected, state=evaluator.state.values)
    try:
        geometry = reconstruct_clipped(evaluator.pipeline.grid, density.projected, evaluator.pipeline.free)
        geometry_report = write_reconstruction(geometry, output / "geometry.json")
        wall_audit = audit_wall_displacement(evaluator.pipeline.grid, density.projected, evaluator.pipeline.free, geometry)
    except (RuntimeError, ValueError, AssertionError) as error:
        geometry_report = None
        wall_audit = {"error": f"{type(error).__name__}: {error}"}
    audit = audit_kkt(design, value.gradient, inequalities, value.jacobian) if stationarity else None
    return {
        "stage": asdict(evaluator.stage), "pde": evaluator.pde, "design": design,
        "raw_density": evaluator.pipeline.raw(design),
        "filtered_density": density.filtered, "projected_density": density.projected, "field_archive": "fields.npz",
        "objective": objective, "power_w_per_m": power, "side_fraction": fraction, "flux": flux, "inequalities": inequalities,
        "projected_gray_measure": float(4. * np.mean(density.projected * (1. - density.projected))),
        "mass_imbalance": abs(flux[1] + flux[2] - flux[0]) / flux[0],
        "nonlinear_residual": value.nonlinear_residual,
        "split_feasible": bool(max(inequalities[:2]) <= FEASIBILITY_TOLERANCE), "kkt": None if audit is None else asdict(audit),
        "geometry": geometry_report, "wall_geometry_audit": wall_audit,
    }


def phase_one(evaluator, incoming, work):
    initial = work.measure("incoming_state", evaluator, lambda: evaluator.evaluate(incoming))
    if max(initial.constraints) <= FEASIBILITY_TOLERANCE:
        return incoming, initial, []
    attempts = []
    for margin in (.05, 0.):
        result = work.measure("slsqp_restoration", evaluator, lambda: restore(
            evaluator, incoming, margin=margin, max_evaluations=RESTORE_BUDGET))
        attempts.append(asdict(result))
        if result.accepted_design is not None:
            return np.asarray(result.accepted_design), initial, attempts
    return None, initial, attempts


def in_box(design, size):
    return (design is not None and design.shape == (size,) and np.isfinite(design).all()
            and np.all((0. <= design) & (design <= 1.)))


def run(output, common, pde="navier-stokes"):
    reference, reference_record = sharp_reference(common)
    grid = DesignGrid(.25)
    domain = build_domain(grid, .25)
    domain.fem.inertia = REYNOLDS
    pipeline = WallDisplacementDensityV1(grid)
    original = pipeline.initial.copy()
    if grid.size != 960 or original.size != np.count_nonzero(pipeline.free):
        raise RuntimeError("The anchored model does not match the original 960-cell design grid")
    if not in_box(original, original.size) or not np.array_equal(pipeline.raw(original), grid.base_raw):
        raise RuntimeError("The anchored initialization did not reproduce the original sharp T")
    save(output, "sharp_reference", reference_record)
    save(output, "setup", {
        "slurm_job_id": os.environ["SLURM_JOB_ID"], "design_cells": grid.size, "free_controls": original.size,
        "free_cells": pipeline.free, "design_spacing_mm": grid.spacing_mm, "fem_spacing_mm": domain.max_spacing_mm,
        "fem_triangles": len(domain.fem.triangles), "fem_dofs": domain.fem.size, "reynolds_number": REYNOLDS,
        "pde": pde, "pde_convective_inertia": 0. if pde == "stokes" else REYNOLDS,
        "flow_band": SIDE_BAND, "flow_constraint_scale": SPLIT_SCALE,
        "normalization_W_per_m": reference, "original_raw_sharp_t_design": grid.base_raw,
        "projection_threshold": .5, "stages": [asdict(stage) for stage in STAGES],
        "restoration_margins": (.05, 0.), "restoration_budget_per_attempt": RESTORE_BUDGET,
        "mma_unique_evaluation_budget_per_stage": MMA_BUDGET,
        "filter_only_reference": "Original sharp raw controls passed through the anchored density map at each stage",
        "geometric_checks": "Recorded independently without entering optimization constraints",
    })

    work = WorkLog()
    incoming = original.copy()
    previous_state = filter_state = None
    stage_results = []
    for index, stage in enumerate(STAGES):
        folder = output / f"stage_{index:02d}"
        folder.mkdir()
        filter_folder = folder / "filter_only"
        filter_folder.mkdir()
        fixed = WallDisplacementEvaluatorV1(domain, pipeline, stage, reference, initial_state=filter_state, pde=pde)
        stage_summary = {"stage": asdict(stage), "qualified_checkpoint": None}
        offset = len(work.entries)
        try:
            baseline = work.measure("filter_only_reference", fixed, lambda: record_design(
                fixed, original, filter_folder, stationarity=False))
            save(filter_folder, "checkpoint", baseline)
            filter_state = fixed.state
            stage_summary["filter_only_checkpoint"] = filter_folder.relative_to(output).as_posix()
            seed = filter_state if previous_state is None else previous_state
            evaluator = WallDisplacementEvaluatorV1(domain, pipeline, stage, reference, initial_state=seed, pde=pde)
            start, initial, restoration = phase_one(evaluator, incoming, work)
            save(folder, "restoration", {
                "stage": asdict(stage), "incoming_design": incoming, "incoming_inequalities": initial.constraints,
                "incoming_side_fraction": initial.side_fraction, "attempts": restoration, "accepted_design": start,
            })
            if start is None:
                stage_summary["stop_reason"] = "SLSQP restoration found no feasible start"
            else:
                outcome = work.measure("mma", evaluator, lambda: optimize(
                    evaluator, start, algorithm="mma", max_evaluations=MMA_BUDGET))
                save(folder, "mma", asdict(outcome))
                native = outcome.optimizer_success and outcome.status in (1, 2, 3, 4)
                stage_summary.update(mma_status=outcome.status, mma_message=outcome.message, native_converged=native)
                candidates = (("terminal", outcome.terminal_design), ("best_feasible_trial", outcome.best_feasible_design),
                              ("restored_start", start))
                checked = []
                for label, design in candidates:
                    if not in_box(design, original.size) or any(np.array_equal(design, prior) for prior in checked):
                        continue
                    checked.append(design)
                    candidate = folder / label
                    candidate.mkdir()
                    replay = WallDisplacementEvaluatorV1(domain, pipeline, stage, reference, initial_state=seed, pde=pde)
                    try:
                        report = work.measure(f"independent_replay_{label}", replay, lambda: record_design(
                            replay, design, candidate, stationarity=True))
                    except (RuntimeError, ValueError, AssertionError, np.linalg.LinAlgError) as error:
                        save(candidate, "checkpoint", {"stage": asdict(stage), "design": design,
                                                       "error": f"{type(error).__name__}: {error}"})
                        continue
                    qualified = (label == "terminal" and native and report["split_feasible"] and report["kkt"]["verified"])
                    report.update(native_terminal_converged=native, qualified_optimization_terminal=qualified,
                                  power_ratio_vs_filter_only=report["power_w_per_m"] / baseline["power_w_per_m"])
                    save(candidate, "checkpoint", report)
                    if label == "terminal":
                        stage_summary.update(independent_terminal_kkt=report["kkt"]["verified"],
                                             terminal_filter_only_power_ratio=report["power_ratio_vs_filter_only"])
                    if qualified:
                        incoming, previous_state = design.copy(), replay.state
                        stage_summary["qualified_checkpoint"] = candidate.relative_to(output).as_posix()
                if stage_summary["qualified_checkpoint"] is None:
                    stage_summary["stop_reason"] = ("MMA did not report genuine convergence" if not native else
                                                    "MMA terminal failed the independent feasibility/KKT replay")
        except (RuntimeError, ValueError, AssertionError, np.linalg.LinAlgError) as error:
            stage_summary["stop_reason"] = f"{type(error).__name__}: {error}"
        stage_summary["work"] = work.entries[offset:]
        save(folder, "summary", stage_summary)
        stage_results.append(stage_summary)
        if stage_summary["qualified_checkpoint"] is None:
            break

    completed = len(stage_results) == len(STAGES) and stage_results[-1]["qualified_checkpoint"] is not None
    return {"pde": pde, "stage_results": stage_results, "scheduled_stages": len(STAGES), "completed_schedule": completed,
            "final_candidate": stage_results[-1]["qualified_checkpoint"] if completed else None,
            "stop_reason": "continued through all verified optimization stages" if completed else stage_results[-1]["stop_reason"],
            "numerical_work_totals": work.totals(), "physical_verification_pending": completed}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="New per-job evidence directory")
    parser.add_argument("--common", type=Path, required=True, help="Directory containing the verified sharp fitted baseline")
    parser.add_argument("--pde", choices=PDE_MODES, default="navier-stokes")
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run numerical work as a Sol Slurm job")
    args.output.mkdir(parents=True, exist_ok=False)
    write_report(args.output, "environment", environment)
    if not write_report(args.output, "summary", lambda: run(args.output, args.common, args.pde)):
        raise SystemExit(1)
    summary = json.loads((args.output / "summary.json").read_text(encoding="utf-8"))["results"]
    if not summary["completed_schedule"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
