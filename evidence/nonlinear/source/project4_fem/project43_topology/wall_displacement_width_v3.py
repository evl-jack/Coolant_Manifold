"""Sol experiment for the protected-width wall displacement with flow-only MMA."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from hashlib import sha256
import json
import os
from pathlib import Path
from time import perf_counter
from unittest.mock import patch

import numpy as np
from shapely.geometry import box
from shapely.ops import unary_union

from . import wall_displacement_geometry_v1 as geometry_tools
from .density import DesignGrid
from .domain import build_domain
from .evaluator import REYNOLDS, Stage
from .full_density_continuation import fitted_screen
from .gate2 import write_report
from .optimize import optimize
from .relaxed_flow import save
from .restoration import FEASIBILITY_TOLERANCE, restore
from .verify_density import STATE_TOLERANCE
from .verify_flow import FINE_STEPS, RESIDUAL_DERIVATIVE_LIMIT, _relative_error
from .wall_displacement_fem_v1 import WallDisplacementEvaluatorV1
from .wall_displacement_optimize_v1 import (
    MMA_BUDGET, RESTORE_BUDGET, SIDE_BAND, SPLIT_SCALE, WorkLog, record_design, sharp_reference,
)
from .wall_displacement_smooth_v2 import GeometricConstraint, ProtectedDensity
from .wall_displacement_verify_v1 import ReferenceDensity, ReferenceFlow, verify


STAGE = Stage(.10, 8., 1e5)
SPACING_MM = .25
SEED_JOB = 2406318
COMPARISON_JOB = 2407390


def setup(checkpoint, common):
    payload = checkpoint.read_bytes()
    archive = json.loads(payload)
    record = archive["results"]
    seed_setup_path = checkpoint.parents[2] / "setup.json"
    seed_setup_payload = seed_setup_path.read_bytes()
    seed_setup = json.loads(seed_setup_payload)
    seed_metadata = seed_setup["results"]
    # The original Navier–Stokes archive predates the explicit PDE field.
    if (not archive["completed"] or not seed_setup["completed"]
            or str(seed_metadata["slurm_job_id"]) != str(SEED_JOB)
            or seed_metadata["reynolds_number"] != REYNOLDS
            or record.get("pde", seed_metadata.get("pde", "navier-stokes")) != "navier-stokes"
            or (record["stage"]["q"], record["stage"]["beta"]) != (.10, 4.)):
        raise ValueError("Expected the completed Navier–Stokes q=0.10, beta=4 archive")
    grid = DesignGrid(SPACING_MM)
    pipeline = ProtectedDensity(grid)
    prior = np.asarray(record["raw_density"], dtype=float)
    if prior.shape != (grid.size,) or not np.all((0. <= prior) & (prior <= 1.)):
        raise ValueError("Invalid archived full raw density")
    reference, provenance = sharp_reference(common)
    metadata = {"seed_archive": str(checkpoint), "seed_sha256": sha256(payload).hexdigest(),
                "seed_job": SEED_JOB, "comparison_job": COMPARISON_JOB,
                "seed_setup_archive": str(seed_setup_path), "seed_setup_sha256": sha256(seed_setup_payload).hexdigest(),
                "power_reference": provenance, "normalization_W_per_m": reference,
                "stage": asdict(STAGE), "reynolds": REYNOLDS, "flow_band": SIDE_BAND,
                "design_cell_mm": SPACING_MM, "max_wall_displacement_mm": 1.,
                "filter_radius_mm": pipeline.filter_radius_mm, "core_boxes_mm": pipeline.core_boxes,
                "editable_cell_count": pipeline.size, "editable_cell_ids": np.flatnonzero(pipeline.free)}
    return pipeline, prior[pipeline.free], reference, metadata


def evaluator_for(pipeline, reference, initial_state=None):
    domain = build_domain(pipeline.grid, SPACING_MM)
    domain.fem.inertia = REYNOLDS
    return WallDisplacementEvaluatorV1(domain, pipeline, STAGE, reference, initial_state=initial_state)


def corridor_checks(pipeline, density, geometry):
    """Check material anchors and the reconstructed port-connected corridor."""
    grid = pipeline.grid
    core = unary_union([box(*bounds) for bounds in pipeline.core_boxes])
    x, y = grid.centers_mm.T
    protected = np.zeros(grid.size, dtype=bool)
    half = grid.spacing_mm / 2.
    for xmin, ymin, xmax, ymax in pipeline.core_boxes:
        protected |= ((x + half > xmin) & (x - half < xmax)
                      & (y + half > ymin) & (y - half < ymax))
    fixed = ~pipeline.free
    anchors = grid.base_raw[fixed]
    checks = {
        "protected_cells_fixed_fluid": bool(not np.any(pipeline.free[protected])
                                            and np.all(grid.base_raw[protected] == 1.)),
        "raw_anchors_exact": bool(np.array_equal(density.raw[fixed], anchors)),
        "filtered_anchors_exact": bool(np.array_equal(density.filtered[fixed], anchors)),
        "projected_anchors_exact": bool(np.array_equal(density.projected[fixed], anchors)),
        "resistance_anchors_exact": bool(np.array_equal(density.resistance[fixed], STAGE.alpha_max * (1. - anchors))),
        "editable_footprint_disjoint_from_corridor": bool(
            geometry_tools.editable_geometry(grid, pipeline.free).intersection(core).area == 0.),
        "port_component_covers_corridor": bool(geometry.port_component is not None
                                               and geometry.port_component.covers(core)),
    }
    return {"passed": all(checks.values()), "checks": checks, "protected_grid_cells": int(protected.sum()),
            "main_corridor_width_mm": 3.75, "side_corridor_width_mm": 1.}


def check_gradient_work(evaluator, design):
    """Exercise the production callback with geometry routines forbidden."""
    with (patch.object(geometry_tools, "reconstruct", side_effect=AssertionError("Geometry in gradient callback")),
          patch.object(GeometricConstraint, "gradient", side_effect=AssertionError("Geometric finite difference")),
          patch.object(geometry_tools, "_turning_report", side_effect=AssertionError("Angle in gradient callback"))):
        value = evaluator.evaluate(design)
        cached = evaluator.evaluate(design) is value
    size = evaluator.pipeline.size
    counts = {"forward_attempts": evaluator.forward_attempts, "forward_solves": evaluator.forward_solves,
              "fluid_adjoint_rhs": evaluator.adjoint_rhs, "geometric_fd_evaluations": 0}
    return {"passed": bool(cached and tuple(counts.values()) == (1, 1, 2, 0)
                           and value.gradient.shape == (size,) and value.jacobian.shape == (2, size)
                           and value.constraints.shape == (2,)),
            **counts, "cached_repeat": cached, "objective_gradient_shape": value.gradient.shape,
            "flow_jacobian_shape": value.jacobian.shape}


def check_complete_residual(evaluator, design):
    """Perturb state and compact design together in the complete Re=50 residual."""
    pipeline, state = evaluator.pipeline, evaluator.state
    fem, free = evaluator.domain.fem, state.free
    reference = ReferenceDensity(pipeline.grid, protected_v2=True)
    flow = ReferenceFlow(evaluator.domain, reference, "navier-stokes")
    rng = np.random.default_rng(43579)
    direction = rng.uniform(-.8, .8, pipeline.size)
    state_direction = np.zeros(fem.size)
    state_direction[free] = rng.normal(size=free.size)
    covector = rng.normal(size=free.size)
    _, _, _, projection_slope, resistance_slope = reference.fields(design, STAGE)
    design_action = flow.residual_direction(state, reference.filter_direction(direction)
                                            * projection_slope * resistance_slope)
    tangent = state.tangent[free][:, free]
    action = tangent @ state_direction[free] + design_action
    transpose_error = _relative_error(covector @ (tangent @ state_direction[free]),
                                      (tangent.T @ covector) @ state_direction[free])

    def perturbed(step):
        resistance = evaluator.transfer.forward(evaluator.density(design + step * direction).resistance)
        operator = fem.viscous_pressure + fem.drag_matrix(resistance)
        return fem.residual(state.values + step * state_direction, operator, evaluator.load, tangent=False)[0][free]

    errors = [{"step": h, "relative_error": _relative_error((perturbed(h) - perturbed(-h)) / (2. * h), action)}
              for h in FINE_STEPS]
    return {"passed": bool(transpose_error <= STATE_TOLERANCE and all(
                sample["relative_error"] <= RESIDUAL_DERIVATIVE_LIMIT for sample in errors[-2:])),
            "state_tangent_transpose_error": transpose_error, "finite_differences": errors}


def check_width(pipeline, seed, reference):
    evaluator = evaluator_for(pipeline, reference)
    interior = .14 + .72 * pipeline.initial
    gradient_work = check_gradient_work(evaluator, interior)
    residual = check_complete_residual(evaluator, interior)
    corridors = {}
    for label, design in (("zero_controls", np.zeros(pipeline.size)), ("sharp_seed", pipeline.initial),
                          ("all_one_controls", np.ones(pipeline.size)), ("remapped_checkpoint", seed)):
        density = evaluator.density(design)
        geometry = geometry_tools.reconstruct_clipped(pipeline.grid, density.projected, pipeline.free)
        corridors[label] = corridor_checks(pipeline, density, geometry)
    return gradient_work, residual, corridors


def verify_width(output, pipeline, seed, reference, common):
    with ThreadPoolExecutor(max_workers=2) as pool:
        fluid_task = pool.submit(verify, common / "sharp_fitted_baseline.json", protected_v2=True)
        width_task = pool.submit(check_width, pipeline, seed, reference)
        fluid = fluid_task.result()
        save(output, "independent_fluid_verification", fluid)
        gradient_work, residual, corridors = width_task.result()
    return {"passed": bool(fluid["passed"] and gradient_work["passed"] and residual["passed"]
                           and all(item["passed"] for item in corridors.values())),
            "sol_job_id": os.environ["SLURM_JOB_ID"], "editable_cell_count": pipeline.size,
            "independent_fluid_verification_passed": fluid["passed"], "gradient_work": gradient_work,
            "complete_residual": residual, "corridors": corridors}


def archive_design(output, label, evaluator, design, *, stationarity=False):
    folder = output / label
    folder.mkdir()
    record = record_design(evaluator, design, folder, stationarity=stationarity)
    density = evaluator.density(design)
    pipeline = evaluator.pipeline
    geometry = geometry_tools.reconstruct_clipped(pipeline.grid, density.projected, pipeline.free)
    record["protected_corridor"] = corridor_checks(pipeline, density, geometry)
    record["turning_angle_role"] = "Archived geometry diagnostic"
    save(folder, "checkpoint", record)
    np.save(folder / "brinkman_resistance.npy", density.resistance, allow_pickle=False)
    return record, geometry


def progress_callback(journal, evaluator, phase, *, margin=None):
    start = perf_counter()
    forward_start, adjoint_start = evaluator.forward_seconds, evaluator.adjoint_seconds
    best = None

    def report(point, unique, total):
        nonlocal best
        violation = max(0., max(point.inequalities))
        best = violation if best is None else min(best, violation)
        restoring = margin is not None
        entry = {"phase": phase, "margin": margin, "elapsed_seconds": perf_counter() - start,
                 "unique_evaluations": unique, "total_evaluations": total,
                 "objective": point.physical_objective if restoring else point.objective,
                 "side_fraction": SIDE_BAND[0] - SPLIT_SCALE * point.inequalities[0],
                 "epigraph_slack": point.slack if restoring else None,
                 "required_feasibility_slack": violation,
                 "best_feasibility_violation": best, "feasible_seen": best <= FEASIBILITY_TOLERANCE,
                 "forward_seconds": evaluator.forward_seconds - forward_start,
                 "adjoint_seconds": evaluator.adjoint_seconds - adjoint_start}
        line = json.dumps(entry, allow_nan=False)
        print(line, file=journal, flush=True)
        print(line, flush=True)

    return report


def run_phases(output, evaluator, seed, work, journal):
    attempts = []
    chosen = None
    for index, margin in enumerate((.05, 0.)):
        progress = progress_callback(journal, evaluator, "phase_one", margin=margin)
        result = work.measure("slsqp_restoration", evaluator, lambda: restore(
            evaluator, seed, margin=margin, max_evaluations=RESTORE_BUDGET, progress=progress))
        attempts.append(result)
        save(output, f"restoration_{index}", asdict(result))
        if result.accepted_design is not None:
            chosen = np.asarray(result.accepted_design)
            break
    if chosen is None:
        return None, attempts
    progress = progress_callback(journal, evaluator, "phase_two")
    return work.measure("mma", evaluator, lambda: optimize(
        evaluator, chosen, algorithm="mma", max_evaluations=MMA_BUDGET, progress=progress)), attempts


def run(output, pipeline, seed, reference):
    evaluator = evaluator_for(pipeline, reference)
    work = WorkLog()
    seed_record, _ = work.measure("archive_seed", evaluator, lambda: archive_design(output, "seed", evaluator, seed))
    seed_state = evaluator.state
    with (output / "progress.jsonl").open("x", encoding="utf-8") as journal:
        outcome, attempts = run_phases(output, evaluator, seed, work, journal)
    if outcome is None:
        trials = (point for attempt in attempts for point in attempt.trial_history)
        best = min(trials, key=lambda point: (point.engineering_violation, point.physical_objective), default=None)
        side = seed_record["side_fraction"] if best is None else SIDE_BAND[0] - SPLIT_SCALE * best.inequalities[0]
        return {"passed": False, "stop_reason": "Phase I flow restoration failed", "work": work.totals(),
                "best_achieved": {"side_fraction": side, "main_fraction": 1. - side,
                                  "trial": None if best is None else asdict(best)}}
    save(output, "mma", asdict(outcome))
    terminal = np.asarray(outcome.terminal_design)
    replay = evaluator_for(pipeline, reference, initial_state=seed_state)
    record, geometry = work.measure("archive_terminal", replay, lambda: archive_design(
        output, "terminal", replay, terminal, stationarity=True))
    corridor = record["protected_corridor"]["passed"]
    fitted = fitted_screen(geometry, SIDE_BAND, area_cap_mm2=None) if corridor else {}
    independently_stationary = record["kkt"]["verified"]
    accepted = (outcome.optimizer_success or outcome.status == 5) and independently_stationary
    return {"passed": bool(accepted and record["split_feasible"] and corridor),
            "mma_status": outcome.status, "mma_message": outcome.message,
            "mma_native_converged": outcome.optimizer_success, "independently_stationary": independently_stationary,
            "objective": record["objective"], "side_fraction": record["side_fraction"],
            "protected_corridor": record["protected_corridor"], "fitted": fitted,
            "work": work.totals(), "work_by_operation": work.entries}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("verify", "optimize"))
    parser.add_argument("output", type=Path)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--common", type=Path, required=True)
    parser.add_argument("--verification", type=Path)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run numerical work in a Sol Slurm job")
    pipeline, seed, reference, metadata = setup(args.checkpoint, args.common)
    if args.mode == "optimize":
        if args.verification is None:
            parser.error("Optimization requires the passing verification summary")
        verification = json.loads(args.verification.read_bytes())
        if not verification["completed"] or not verification["results"]["passed"]:
            parser.error("Width-only verification did not pass")
        verified_setup = json.loads((args.verification.parent / "setup.json").read_bytes())["results"]
        if (verified_setup["seed_sha256"] != metadata["seed_sha256"]
                or verified_setup["seed_setup_sha256"] != metadata["seed_setup_sha256"]
                or verified_setup["power_reference"]["sha256"] != metadata["power_reference"]["sha256"]
                or verified_setup["editable_cell_ids"] != metadata["editable_cell_ids"].tolist()):
            parser.error("Verification inputs differ from optimization inputs")
        metadata["verification_summary"] = str(args.verification)
    args.output.mkdir(parents=True, exist_ok=False)
    save(args.output, "setup", metadata)
    if args.mode == "verify":
        task = lambda: verify_width(args.output, pipeline, seed, reference, args.common)
    else:
        task = lambda: run(args.output, pipeline, seed, reference)
    completed = write_report(args.output, "summary", task)
    passed = completed and json.loads((args.output / "summary.json").read_bytes())["results"]["passed"]
    raise SystemExit(0 if passed else 2)


if __name__ == "__main__":
    main()
