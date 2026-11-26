"""Screen reduced-basis flow authority and fit an archived full-space anchor on Sol."""

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
from time import perf_counter
from types import SimpleNamespace

import numpy as np
from scipy.optimize import lsq_linear

from .density import DensityPipeline, DesignGrid, ReducedDensity
from .domain import build_domain
from .evaluator import Evaluator, Stage, baseline_reference
from .kkt_audit import audit_kkt
from .optimize import optimize


STAGE = Stage(0.01, 0., 1e5)
DESIGN_SPACING_MM = 0.5
ANALYSIS_SPACING_MM = 0.5
FEASIBILITY_TOLERANCE = 1e-6
START_SEED = 43205


class AreaOnlyFlow:
    """Minimize minus side-flow fraction with only the area inequality."""

    def __init__(self, evaluator):
        self.evaluator = evaluator

    def evaluate(self, design):
        evaluation = self.evaluator.evaluate(design)
        return SimpleNamespace(objective=-evaluation.side_fraction,
                               gradient=-0.001 * evaluation.jacobian[1],
                               constraints=evaluation.constraints[2:],
                               jacobian=evaluation.jacobian[2:])


def _starts():
    expansion, contraction = np.r_[np.ones(10), np.zeros(10)], np.r_[np.zeros(10), np.ones(10)]
    archived_corner = np.array([1., 1., 0., 0., 0., 1., 1., 0., 0., 0.] + [1.] * 10)
    rng = np.random.default_rng(START_SEED)
    return [("zero", np.zeros(20)), ("all_one", np.ones(20)), ("expansion_only", expansion),
            ("contraction_only", contraction), ("archived_corner", archived_corner),
            *((f"seeded_interior_{index}", rng.uniform(0.05, 0.95, 20)) for index in range(1, 4))]


def _responses(evaluation):
    return {"eta": evaluation.side_fraction, "eta_pct": 100 * evaluation.side_fraction,
            "area_mm2": evaluation.area_mm2, "power_w_per_m": evaluation.power_w_per_m,
            "normalized_power": evaluation.objective, "engineering_inequalities": evaluation.constraints.tolist(),
            "area_feasible": bool(evaluation.constraints[2] <= FEASIBILITY_TOLERANCE),
            "engineering_feasible": bool(max(evaluation.constraints) <= FEASIBILITY_TOLERANCE),
            "mass_imbalance": evaluation.mass_imbalance, "nonlinear_residual": evaluation.nonlinear_residual,
            "gray_measure": evaluation.gray_measure, "flux": evaluation.flux.tolist()}


def _cost(evaluator):
    return {counter: getattr(evaluator, counter) for counter in
            ("forward_attempts", "forward_solves", "adjoint_rhs", "newton_iterations")}


def _replay_endpoint(make_evaluator, design):
    evaluator = make_evaluator()
    try:
        value = evaluator.evaluate(design)
        area_objective = -0.001 * value.jacobian[1]
        witness = asdict(audit_kkt(design, area_objective, value.constraints[2:], value.jacobian[2:]))
        witness.update(multipliers=witness["multipliers"].tolist(),
                       projected_residual=witness["projected_residual"].tolist())
        return {"design": design.tolist(), "response": _responses(value), "cost": _cost(evaluator),
                "kkt": witness, "objective_gradient": area_objective.tolist(),
                "area_jacobian": value.jacobian[2].tolist()}
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
        return {"design": design.tolist(), "error": f"{type(error).__name__}: {error}", "cost": _cost(evaluator)}


def _screen_start(name, initial, make_evaluator, max_evaluations):
    start = perf_counter()
    evaluator = make_evaluator()
    outcome = optimize(AreaOnlyFlow(evaluator), initial, algorithm="SLSQP", max_evaluations=max_evaluations)
    terminal = _replay_endpoint(make_evaluator, outcome.terminal_design)
    best = None
    if outcome.best_feasible_design is not None:
        if np.array_equal(outcome.best_feasible_design, outcome.terminal_design):
            best = terminal
        else:
            best = _replay_endpoint(make_evaluator, outcome.best_feasible_design)
    history = [{"design": item.design.tolist(), "eta": -item.objective,
                "area_mm2": 83. + float(item.inequalities[0]), "area_feasible": item.feasible}
               for item in outcome.evaluation_history]
    return {"name": name, "initial_design": initial.tolist(), "terminal": terminal,
            "terminal_source": outcome.terminal_source, "best_observed_area_feasible": best,
            "optimizer_success": outcome.optimizer_success, "status": outcome.status, "message": outcome.message,
            "iterations": outcome.iterations, "objective_calls": outcome.objective_calls,
            "unique_evaluations": outcome.unique_evaluations, "optimizer_cost": _cost(evaluator),
            "history": history, "wall_seconds": perf_counter() - start}


def _coverage(grid, pipeline, reduced):
    base = grid.base_raw
    # These extrema are simultaneous. All contraction columns are nonpositive
    # and all expansion columns nonnegative, so one corner attains each vector.
    raw_minimum = base + np.minimum(reduced.jacobian, 0.).sum(axis=1)
    raw_maximum = base + np.maximum(reduced.jacobian, 0.).sum(axis=1)
    minimum = pipeline.evaluate(raw_minimum, **asdict(STAGE))
    maximum = pipeline.evaluate(raw_maximum, **asdict(STAGE))
    baseline = pipeline.evaluate(base, **asdict(STAGE))
    expansion = raw_maximum > base
    contraction = raw_minimum < base
    return {"shape_ny_nx": list(grid.shape), "x_edges_mm": grid.x_edges_mm.tolist(),
            "y_edges_mm": grid.y_edges_mm.tolist(), "ordering": "y outer, x inner",
            "scope": "Pointwise box-only extrema. Does not impose the area or flow constraints.",
            "raw_baseline": base.tolist(), "raw_minimum": raw_minimum.tolist(), "raw_maximum": raw_maximum.tolist(),
            "filtered_baseline": baseline.filtered.tolist(), "filtered_minimum": minimum.filtered.tolist(),
            "filtered_maximum": maximum.filtered.tolist(), "projected_baseline": baseline.projected.tolist(),
            "projected_minimum": minimum.projected.tolist(), "projected_maximum": maximum.projected.tolist(),
            "raw_expansion_coverage_mask": expansion.tolist(),
            "raw_contraction_coverage_mask": contraction.tolist(),
            "covered_expansion_cells": int(np.count_nonzero(expansion)),
            "covered_contraction_cells": int(np.count_nonzero(contraction)),
            "raw_expansion_maximum_change": (raw_maximum - base).tolist(),
            "raw_contraction_maximum_change": (base - raw_minimum).tolist(),
            "filtered_minimum_area_mm2": minimum.area_mm2, "filtered_maximum_area_mm2": maximum.area_mm2}


def _fit_anchor(anchor_path, make_full, make_reduced, pipeline, reduced):
    report = json.loads(anchor_path.read_text(encoding="utf-8"))
    summary = report.get("results", report)
    if not summary or not summary.get("restored") or "anchor" not in summary:
        raise ValueError("The supplied restoration summary has no accepted full-space anchor")
    raw = np.asarray(summary["anchor"], dtype=float)
    if raw.shape != (pipeline.grid.size,) or not np.all(np.isfinite(raw)) or np.any((raw < 0) | (raw > 1)):
        raise ValueError("Expected a length-240 unit-box full-space anchor on the common grid")
    fit = lsq_linear(reduced.jacobian, raw - pipeline.grid.base_raw, bounds=(0.0, 1.0), tol=1e-12)
    original_evaluator, fitted_evaluator = make_full(), make_reduced()
    original, fitted = original_evaluator.evaluate(raw), fitted_evaluator.evaluate(fit.x)
    original_density, fitted_density = original_evaluator.density(raw), fitted_evaluator.density(fit.x)
    errors = {}
    for field in ("raw", "filtered", "projected"):
        difference = getattr(fitted_density, field) - getattr(original_density, field)
        errors[field] = {"rmse": float(np.sqrt(np.mean(difference**2))),
                         "linf": float(np.max(np.abs(difference)))}
    return {"source": str(anchor_path), "initial_archived_inequalities": summary.get("anchor_constraints"),
            "fit_success": bool(fit.success), "fit_status": int(fit.status), "fit_message": str(fit.message),
            "fit_optimality": float(fit.optimality), "fit_cost": float(fit.cost), "fitted_controls": fit.x.tolist(),
            "density_error": errors, "original_raw": raw.tolist(), "fitted_raw": fitted_density.raw.tolist(),
            "original_projected": original_density.projected.tolist(), "fitted_projected": fitted_density.projected.tolist(),
            "original_response": _responses(original), "fitted_response": _responses(fitted),
            "original_cost": _cost(original_evaluator), "fitted_cost": _cost(fitted_evaluator),
            "delta_eta_pp_fitted_minus_original": 100 * (fitted.side_fraction - original.side_fraction),
            "delta_area_mm2_fitted_minus_original": fitted.area_mm2 - original.area_mm2,
            "delta_power_w_per_m_fitted_minus_original": fitted.power_w_per_m - original.power_w_per_m}


def _plot(output, coverage, fit):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fields = [("Original sharp T-junction", coverage["raw_baseline"]),
              ("Minimum attainable raw density", coverage["raw_minimum"]),
              ("Maximum attainable raw density", coverage["raw_maximum"]),
              ("Minimum attainable projected density", coverage["projected_minimum"]),
              ("Maximum attainable projected density", coverage["projected_maximum"])]
    if "error" not in fit:
        fields.extend([("Archived full-space anchor", fit["original_raw"]),
                       ("Reduced fit to anchor", fit["fitted_raw"]),
                       ("Absolute raw density error", np.abs(np.subtract(fit["original_raw"], fit["fitted_raw"])))])
    fig, axes = plt.subplots(2, 4, figsize=(13, 6), layout="constrained")
    for axis, (title, cells) in zip(axes.ravel(), fields):
        image = axis.pcolormesh(coverage["x_edges_mm"], coverage["y_edges_mm"],
                                np.reshape(cells, coverage["shape_ny_nx"]), vmin=0, vmax=1, cmap="viridis")
        axis.set(title=title, xlabel="x (mm)", ylabel="y (mm)", aspect="equal")
    for axis in axes.ravel()[len(fields):]:
        axis.set_visible(False)
    fig.colorbar(image, ax=axes.ravel().tolist(), label="Density")
    fig.savefig(output / "reduced_authority_maps.png", dpi=180)
    plt.close(fig)


def _write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def run(output, full_anchor, max_evaluations, plot=False):
    grid = DesignGrid(DESIGN_SPACING_MM)
    pipeline = DensityPipeline(grid)
    reduced = ReducedDensity(grid)
    domain = build_domain(grid, ANALYSIS_SPACING_MM)
    p0, seed = baseline_reference(domain, pipeline, STAGE.alpha_max)

    def full_evaluator():
        return Evaluator(domain, pipeline, STAGE, p0, initial_state=seed)

    def reduced_evaluator():
        return Evaluator(domain, pipeline, STAGE, p0, reduced=reduced, initial_state=seed)

    coverage = _coverage(grid, pipeline, reduced)
    _write_json(output / "density_coverage.json", coverage)
    runs = []
    for name, initial in _starts():
        result = _screen_start(name, initial, reduced_evaluator, max_evaluations)
        _write_json(output / f"{name}.json", result)
        print(f"Finished {name}: {result['message']}", flush=True)
        runs.append({key: value for key, value in result.items() if key != "history"})
    try:
        fit = _fit_anchor(full_anchor, full_evaluator, reduced_evaluator, pipeline, reduced)
    except (OSError, ValueError, RuntimeError, np.linalg.LinAlgError) as error:
        fit = {"source": str(full_anchor), "error": f"{type(error).__name__}: {error}"}
    _write_json(output / "full_anchor_fit.json", fit)
    plot_error = None
    if plot:
        try:
            _plot(output, coverage, fit)
        except ImportError as error:
            plot_error = str(error)
    observations = [(item["name"], item["best_observed_area_feasible"]) for item in runs
                    if item["best_observed_area_feasible"] and "response" in item["best_observed_area_feasible"]
                    and item["best_observed_area_feasible"]["response"]["area_feasible"]]
    best_observation = max(observations, key=lambda item: item[1]["response"]["eta"], default=None)
    stationary = [item["name"] for item in runs if item["terminal"].get("kkt", {}).get("verified", False)]
    return {"stage": asdict(STAGE), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            "analysis_spacing_mm": ANALYSIS_SPACING_MM, "design_spacing_mm": DESIGN_SPACING_MM,
            "design_cells": grid.size, "triangles": len(domain.fem.triangles), "dofs": domain.fem.size,
            "baseline_power_w_per_m": p0, "start_seed": START_SEED,
            "max_evaluations_per_start": max_evaluations, "starts": runs,
            "best_observed_area_feasible_run_and_replay": best_observation,
            "kkt_verified_terminal_names": stationary, "fit_report": "full_anchor_fit.json",
            "fit_solved": fit.get("fit_success", False),
            "original_anchor_replay_feasible": fit.get("original_response", {}).get("engineering_feasible", False),
            "plot_error": plot_error}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("full_anchor", type=Path, help="Path to full/restoration_summary.json")
    parser.add_argument("--max-evaluations", type=int, default=250)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this diagnostic in a Sol Slurm job")
    if args.max_evaluations < 1:
        parser.error("--max-evaluations must be positive")
    args.output.mkdir(parents=True)
    result = run(args.output, args.full_anchor, args.max_evaluations, args.plot)
    _write_json(args.output / "reduced_authority_summary.json", result)
    if (not result["fit_solved"] or not result["original_anchor_replay_feasible"]
            or not result["best_observed_area_feasible_run_and_replay"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
