"""Independently audit archived wall-displacement numerical evidence on Sol."""

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from .density import DesignGrid
from .domain import build_domain
from .evaluator import Stage
from .verify_density import STATE_TOLERANCE
from .verify_flow import MASS_LIMIT, RESIDUAL_LIMIT, _power_and_flux
from .wall_displacement_density_v1 import WallDisplacementDensityV1
from .wall_displacement_fem_v1 import WallDisplacementEvaluatorV1
from .wall_displacement_verify_v1 import ReferenceDensity, ReferenceFlow


PRIMAL = 1e-6
DUAL = 1e-8
STATIONARITY = 1e-4
COMPLEMENTARITY = 1e-5
MMA_BUDGET = 500
RESTORE_BUDGET = 250
SCHEDULE = ((.01, 0.), (.05, 0.), (.10, 0.), (.10, 1.), (.10, 2.), (.10, 4.), (.10, 8.))


def read(path):
    data = json.loads(path.read_bytes())
    if not data.get("completed"):
        raise ValueError(f"Archive {path} did not complete")
    return data["results"]


def _same(a, b, atol=1e-12):
    return bool(np.allclose(a, b, rtol=1e-10, atol=atol))


def _split(fraction):
    return np.array([(0.10 - fraction) / .001, (fraction - 0.90) / .001])


def _sharp_reference(lane, setup, path):
    reference = read(lane / "sharp_reference.json")
    source = Path(reference["archive"]) if path is None else path
    payload = source.read_bytes()
    fitted = read(source)
    fine = fitted["levels"][-1]
    power = fine["mechanical_power_W_per_m"]
    return {"source": str(source), "sha256": sha256(payload).hexdigest(), "fine": fine,
            "verified": bool(fitted["fitted_diagnostic_passed"] and power > 0
                             and power == setup["normalization_W_per_m"]
                             and fine == reference["finest_level"]
                             and sha256(payload).hexdigest() == reference["sha256"])}


def _witness(evaluation, design, inequality, archived):
    multipliers = np.asarray(archived["multipliers"])
    lagrangian = evaluation.gradient + multipliers @ evaluation.jacobian
    projected = np.minimum(np.maximum(lagrangian, design - 1.), design)
    scalars = dict(primal_violation=float(max(0., inequality.max(), -design.min(), design.max() - 1)),
                   dual_violation=float(max(0., (-multipliers).max())),
                   stationarity_inf=float(np.linalg.norm(projected, np.inf)),
                   complementarity_inf=float(np.max(np.abs(multipliers * inequality))))
    limits = (PRIMAL, DUAL, STATIONARITY, COMPLEMENTARITY)
    verified = all(value <= limit for value, limit in zip(scalars.values(), limits))
    # Separate residual agreement from the threshold classification so a rounding-level
    # disagreement at a limit cannot promote an endpoint to a verified checkpoint.
    agreement = all(_same(value, archived[name], atol=STATE_TOLERANCE) for name, value in scalars.items())
    agreement &= _same(projected, archived["projected_residual"], atol=STATE_TOLERANCE)
    return {**scalars, "multipliers": multipliers.tolist(), "verified": verified,
            "agrees_with_archive": bool(agreement and verified == archived["verified"])}


def _checkpoint(folder, stage, domain, pipeline, reference, flow, normalization, *, baseline=False):
    archived = read(folder / "checkpoint.json")
    if "error" in archived:
        return {"archive_error": archived["error"], "passed": False}
    design = np.asarray(archived["design"])
    expected_raw = reference.raw(design)
    filtered, projected, resistance, _, _ = reference.fields(design, stage)
    with np.load(folder / archived["field_archive"], allow_pickle=False) as arrays:
        state = arrays["state"]
        raw_matches = np.array_equal(arrays["raw"], expected_raw)
        fields_match = (raw_matches and np.array_equal(arrays["compact_design"], design)
                        and _same(arrays["filtered"], filtered, atol=STATE_TOLERANCE)
                        and _same(arrays["projected"], projected, atol=STATE_TOLERANCE))
        fixed_exact = all(np.array_equal(arrays[name][~reference.editable], reference.sharp[~reference.editable])
                          for name in ("raw", "filtered", "projected"))
    fields_match &= np.array_equal(archived["raw_density"], expected_raw)
    fields_match &= _same(archived["filtered_density"], filtered, atol=STATE_TOLERANCE)
    fields_match &= _same(archived["projected_density"], projected, atol=STATE_TOLERANCE)
    fixed_exact &= all(np.array_equal(np.asarray(archived[name])[~reference.editable], reference.sharp[~reference.editable])
                       for name in ("raw_density", "filtered_density", "projected_density"))
    fem = flow.fem
    operator = fem.viscous_pressure + fem.drag_matrix(flow.resistance(resistance))
    residual, _ = fem.residual(state, operator, flow.load, tangent=False)
    nonlinear = float(np.linalg.norm(residual[flow.free], np.inf))
    prescribed = np.array_equal(state[domain.fixed], domain.prescribed)
    power, flux = _power_and_flux(domain, state)
    if not np.isfinite(flux).all() or flux[0] <= 0:
        raise ValueError(f"Nonpositive or nonfinite inlet flux in {folder}")
    fraction = float(flux[2] / flux[0])
    mass = float(abs(flux[0] - flux[1] - flux[2]) / flux[0])
    inequality = _split(fraction)
    feasible = bool(max(inequality) <= PRIMAL)
    response_agrees = (_same(power / normalization, archived["objective"])
                       and _same(inequality, archived["inequalities"], atol=1e-9)
                       and _same(flux, archived["flux"]) and feasible == archived["split_feasible"])
    kkt = None
    if not baseline:
        seed = SimpleNamespace(values=state)
        evaluator = WallDisplacementEvaluatorV1(domain, pipeline, stage, normalization,
                                                  initial_state=seed, pde=flow.pde)
        kkt = _witness(evaluator.evaluate(design), design, inequality, archived["kkt"])
    init_ok = not baseline or np.array_equal(expected_raw, reference.sharp)
    passed = (fields_match and fixed_exact and prescribed and nonlinear <= RESIDUAL_LIMIT
              and mass <= MASS_LIMIT and abs(flux[0] - 1.) <= MASS_LIMIT and response_agrees and init_ok
              and (kkt is None or kkt["agrees_with_archive"]))
    return {"power_w_per_m": float(power), "objective": float(power / normalization), "side_share_pct": 100 * fraction,
            "split_margin_pp": float(100 * min(fraction - 0.10, 0.90 - fraction)), "g": inequality.tolist(),
            "split_feasible": feasible, "mass_imbalance": mass, "nonlinear_residual": nonlinear,
            "fixed_anchors_exact": bool(fixed_exact), "prescribed_velocity_exact": bool(prescribed),
            "fields_agree": bool(fields_match), "responses_agree": bool(response_agrees),
            "sharp_raw_initialization": bool(init_ok) if baseline else None, "kkt": kkt, "passed": bool(passed)}


def _mma(folder):
    outcome = read(folder / "mma.json")
    history = outcome["evaluation_history"]
    unique = len({tuple(point["design"]) for point in history})
    feasibility = all(point["feasible"] == (max(point["inequalities"]) <= PRIMAL) for point in history)
    calls = outcome["objective_calls"]
    native = outcome["status"] in (1, 2, 3, 4) and outcome["optimizer_success"]
    terminal_path = folder / "terminal" / "checkpoint.json"
    same_terminal = (not terminal_path.exists()
                     or outcome["terminal_design"] == read(terminal_path)["design"])
    return {"status": outcome["status"], "native_converged": native, "native_objective_calls": calls,
            "unique_successful_evaluations": unique, "reported_unique_evaluations": outcome["unique_evaluations"],
            "history_feasibility_consistent": feasibility, "terminal_design_agrees": same_terminal,
            "passed": bool(unique == len(history) == outcome["unique_evaluations"] <= MMA_BUDGET
                           and 0 <= calls <= MMA_BUDGET and feasibility and same_terminal)}


def _restoration(folder):
    record = read(folder / "restoration.json")
    attempts = record["attempts"]
    margins = [attempt["margin"] for attempt in attempts]
    budgets = all(attempt["max_evaluations"] == RESTORE_BUDGET
                  and attempt["unique_evaluations"] <= attempt["physical_attempts"] <= RESTORE_BUDGET
                  for attempt in attempts)
    initial = np.asarray(record["incoming_inequalities"])
    initial_agrees = _same(initial, _split(record["incoming_side_fraction"]), atol=1e-9)
    accepted = record["accepted_design"]
    if attempts:
        witnesses = [attempt["accepted"] for attempt in attempts if attempt["accepted"] is not None]
        accepted_valid = (accepted is None and not witnesses) or (
            len(witnesses) == 1 and witnesses[0]["design"] == accepted
            and max(witnesses[0]["inequalities"]) + margins[-1] <= PRIMAL)
    else:
        accepted_valid = accepted == record["incoming_design"] and max(initial) <= PRIMAL
    return {"margins_attempted": margins, "feasibility_accepted": accepted is not None,
            "physical_attempts": [attempt["physical_attempts"] for attempt in attempts],
            "passed": bool(margins in ([], [.05], [.05, 0.]) and budgets and initial_agrees and accepted_valid)}


def _snapshot(lane, live):
    try:
        return read(lane / "summary.json"), True
    except (FileNotFoundError, json.JSONDecodeError):
        if not live:
            raise
    finished = []
    for index in range(len(SCHEDULE)):
        try:
            finished.append(read(lane / f"stage_{index:02d}" / "summary.json"))
        except (FileNotFoundError, json.JSONDecodeError):
            break
    return {"stage_results": finished, "stop_reason": None}, False


def run(lane, sharp_path, *, live=False):
    setup = read(lane / "setup.json")
    summary, final_archive = _snapshot(lane, live)
    mode = setup.get("pde", "navier-stokes")
    grid = DesignGrid(.25)
    domain = build_domain(grid, .25)
    domain.fem.inertia = 0. if mode == "stokes" else 50.
    pipeline = WallDisplacementDensityV1(grid)
    reference = ReferenceDensity(grid)
    flow = ReferenceFlow(domain, reference, mode)
    normalization = setup["normalization_W_per_m"]
    sharp = _sharp_reference(lane, setup, sharp_path)
    setup_ok = (np.array_equal(reference.editable, setup["free_cells"])
                and np.array_equal(reference.sharp, setup["original_raw_sharp_t_design"])
                and setup["free_controls"] == np.count_nonzero(reference.editable)
                and setup["design_cells"] == grid.size
                and setup["flow_band"] == [.10, .90] and setup["flow_constraint_scale"] == .001
                and setup["mma_unique_evaluation_budget_per_stage"] == MMA_BUDGET
                and setup["restoration_budget_per_attempt"] == RESTORE_BUDGET
                and setup["stages"] == [dict(q=q, beta=beta, alpha_max=1e5) for q, beta in SCHEDULE])
    rows = []
    incoming = reference.sharp[reference.editable].tolist()
    for index, result in enumerate(summary["stage_results"]):
        folder = lane / f"stage_{index:02d}"
        stage = Stage(**result["stage"])
        baseline = _checkpoint(folder / "filter_only", stage, domain, pipeline, reference, flow, normalization,
                               baseline=True) if (folder / "filter_only" / "checkpoint.json").exists() else None
        mma = _mma(folder) if (folder / "mma.json").exists() else None
        restoration = _restoration(folder) if (folder / "restoration.json").exists() else None
        candidates = {label: _checkpoint(folder / label, stage, domain, pipeline, reference, flow, normalization)
                      for label in ("terminal", "best_feasible_trial", "restored_start")
                      if (folder / label / "checkpoint.json").exists()}
        terminal = candidates.get("terminal") or {}
        qualified = bool(mma and mma["native_converged"] and terminal.get("split_feasible")
                         and terminal["kkt"]["verified"])
        expected = f"stage_{index:02d}/terminal" if qualified else None
        stage_ok = (result["qualified_checkpoint"] == expected and result == read(folder / "summary.json")
                    and result["stage"] == setup["stages"][index])
        if restoration is not None:
            stage_ok &= read(folder / "restoration.json")["incoming_design"] == incoming
        if expected is not None:
            incoming = read(lane / expected / "checkpoint.json")["design"]
        ratio = None
        if baseline and terminal and "power_w_per_m" in terminal:
            ratio = terminal["power_w_per_m"] / baseline["power_w_per_m"]
            stage_ok &= _same(ratio, result["terminal_filter_only_power_ratio"])
        passed = (stage_ok and baseline is not None and baseline["passed"] and restoration is not None
                  and restoration["passed"] and (mma is None or mma["passed"])
                  and all(item["passed"] for item in candidates.values()))
        rows.append({"index": index, "stage": result["stage"], "filter_only": baseline, "mma": mma,
                     "restoration": restoration, "candidates": candidates, "qualified": qualified,
                     "terminal_to_filter_only_power_ratio": ratio, "passed": bool(passed)})
    schedule_completed = len(rows) == len(SCHEDULE) and all(row["qualified"] for row in rows)
    final_checkpoint = f"stage_{len(rows) - 1:02d}/terminal" if schedule_completed else None
    progression = all(row["qualified"] for row in rows[:-1])
    if final_archive:
        progression &= (bool(rows) and summary["completed_schedule"] == schedule_completed
                        and summary["final_candidate"] == final_checkpoint)
    return {"archive": str(lane), "pde": mode, "source_slurm_job_id": setup["slurm_job_id"],
            "sharp_fitted_reference": sharp, "setup_verified": bool(setup_ok), "stages": rows,
            "lane_summary_available": final_archive, "audited_stage_count": len(rows),
            "continuation_completed": schedule_completed if final_archive else None, "progression_verified": progression,
            "archived_stop_reason": summary["stop_reason"],
            "audit_passed": (bool(setup_ok and sharp["verified"] and progression and all(row["passed"] for row in rows))
                             if rows else None)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lane", type=Path, help="Archived wall-displacement optimization directory")
    parser.add_argument("--sharp-reference", type=Path, help="Override archived path to the same hashed fitted baseline")
    parser.add_argument("--live-snapshot", action="store_true", help="Audit the complete stage prefix of an active run")
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this audit as a Sol Slurm job")
    report = run(args.lane, args.sharp_reference, live=args.live_snapshot)
    print(json.dumps(report, indent=2, allow_nan=False))
    raise SystemExit(1 if report["audit_passed"] is False else 0)


if __name__ == "__main__":
    main()
