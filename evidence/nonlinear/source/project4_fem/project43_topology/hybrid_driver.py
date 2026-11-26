"""Run matched full-density plus explicit-baffle MMA diagnostics on ASU Sol."""

import argparse
from dataclasses import asdict, replace
import json
import os
from pathlib import Path

import numpy as np

from .connectivity import ConnectivityGraph, RectangularBaffleClearance
from .density import DesignGrid
from .domain import build_domain
from .evaluator import Stage
from .full_density_continuation import FITTED_SPACINGS_MM, MAX_AREA_MM2, ROUTES, fitted_screen
from .gate2 import environment, write_report
from .gate2_restoration import FEASIBILITY_TOLERANCE, measured
from .hybrid_density import BAFFLE_X_MM, HEIGHT_LIMIT, HEIGHT_MM, HybridDensityModel, HybridEvaluator, hybrid_reference
from .hybrid_postprocess import export_endpoint, reconstruct_hybrid
from .kkt_audit import audit_kkt
from .optimize import optimize
from .reconstruction import Reconstruction, solve_fitted, write_reconstruction
from .restoration import restore
from .verify_flow import _power_and_flux


BETA_ZERO_STAGE = Stage(0.01, 0., 1e5)
HEIGHT_COSTS = {"model-a": 0., "model-b": 0.1}
BAFFLE_HEIGHTS_MM = (0., .5, 1., 1.5, 2., 2.5, 3., 3.5, 3.75)
MIN_SEED_HEIGHT_MM = 1e-8


def physical_design(unit_design):
    design = np.asarray(unit_design, dtype=float).copy()
    design[-1] *= HEIGHT_LIMIT
    return design


def unit_design(physical):
    design = np.asarray(physical, dtype=float).copy()
    design[-1] /= HEIGHT_LIMIT
    return design


class UnitBoxEvaluator:
    """Pass a unit box to the verified optimizer and apply the exact chain rule in s."""

    def __init__(self, hybrid):
        self.hybrid = hybrid

    def evaluate(self, design):
        response = self.hybrid.evaluate(physical_design(design))
        gradient, jacobian = response.gradient.copy(), response.jacobian.copy()
        gradient[-1] *= HEIGHT_LIMIT
        jacobian[:, -1] *= HEIGHT_LIMIT
        return replace(response, gradient=gradient, jacobian=jacobian)


def preserve(output, name, data):
    if not write_report(output, name, lambda: data):
        raise RuntimeError(f"Could not preserve {name} evidence in {output}")


def response_summary(response):
    return {"power_W_per_m": response.power_w_per_m, "objective": response.objective,
            "side_fraction": response.side_fraction, "area_mm2": response.area_mm2,
            "inequalities": response.constraints, "fluxes": response.flux,
            "mass_imbalance": response.mass_imbalance, "nonlinear_residual": response.nonlinear_residual}


def fitted_baffle_only(height_mm):
    """Use the exact sharp T for the baffle-only sweep, before density filtering."""
    from shapely.geometry import box

    polygon = box(0., 0., 20., 4.).union(box(9.5, 4., 10.5, 7.))
    if height_mm > 0.:
        polygon = polygon.difference(box(BAFFLE_X_MM[0], 0., BAFFLE_X_MM[1], height_mm))
    if polygon.geom_type != "Polygon" or not polygon.is_valid:
        raise RuntimeError("The baffle-only sharp-T face did not form one valid component")
    np.testing.assert_allclose(polygon.area, MAX_AREA_MM2 - (BAFFLE_X_MM[1] - BAFFLE_X_MM[0]) * height_mm,
                               rtol=1e-12)
    geometry = Reconstruction(polygon, (polygon,), (0, 0, 0), ())
    (response, _, _), cost = measured(lambda: solve_fitted(geometry, FITTED_SPACINGS_MM[0]))
    return dict(response, cost=cost)


def baffle_sweep(evaluator, model):
    """Sweep common heights and independently refine the two 50% split brackets."""
    records = []

    def record(height_mm, source):
        design = model.initial.copy()
        design[-1] = height_mm / HEIGHT_MM
        row = dict(height_mm=height_mm, source=source, physical_design=design)
        try:
            value, cost = measured(lambda: evaluator.evaluate(design), evaluator)
            certificate = evaluator.certificate(design)
            geometry = reconstruct_hybrid(model.grid, evaluator.density(design).projected, height_mm)
            row.update(hybrid_immersed=dict(response_summary(value), cost=cost),
                       strict_graph_paths=certificate.hard_paths,
                       reconstructed_connected=geometry.port_component is not None,
                       reconstructed_area_mm2=geometry.area_mm2)
        except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
            row["hybrid_error"] = f"{type(error).__name__}: {error}"
        try:
            row["fitted_sharp_baffle"] = fitted_baffle_only(height_mm)
        except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
            row["fitted_error"] = f"{type(error).__name__}: {error}"
        records.append(row)
        return row

    for height_mm in BAFFLE_HEIGHTS_MM:
        record(height_mm, "fixed_heights")

    def fitted_fraction(row):
        value = row.get("fitted_sharp_baffle", {})
        return value["side_share"] if value.get("state_passed") else None

    def immersed_fraction(row):
        value = row.get("hybrid_immersed")
        if (value is None or not row.get("reconstructed_connected") or not all(row["strict_graph_paths"])
                or row["reconstructed_area_mm2"] > MAX_AREA_MM2
                or max(value["inequalities"][2:]) > FEASIBILITY_TOLERANCE):
            return None
        return value["side_fraction"]

    coarse = records.copy()
    for source, get_fraction in (("fitted_bisection", fitted_fraction), ("immersed_bisection", immersed_fraction)):
        pairs = [(a, b) for a, b in zip(coarse[:-1], coarse[1:])
                 if get_fraction(a) is not None and get_fraction(b) is not None]
        bracket = next(((a, b) for a, b in pairs if (get_fraction(a) - .5) * (get_fraction(b) - .5) < 0.), None)
        if bracket is None:
            continue
        low, high = bracket
        for _ in range(12):
            middle = record((low["height_mm"] + high["height_mm"]) / 2., source)
            fraction = get_fraction(middle)
            if fraction is None or abs(fraction - .5) <= 1e-5:
                break
            if (get_fraction(low) - .5) * (fraction - .5) > 0.:
                low = middle
            else:
                high = middle
    return records


def exact_feasible(hybrid, normalized_design, *, positive=False):
    design = physical_design(normalized_design)
    height_mm = HEIGHT_MM * design[-1]
    if positive and height_mm <= MIN_SEED_HEIGHT_MM:
        return False
    value = hybrid.evaluate(design)
    if max(value.constraints) > FEASIBILITY_TOLERANCE or not all(hybrid.certificate(design).hard_paths):
        return False
    try:
        geometry = reconstruct_hybrid(hybrid.pipeline.grid, hybrid.density(design).projected, height_mm)
    except ValueError:
        return False
    return geometry.port_component is not None and geometry.area_mm2 <= MAX_AREA_MM2


def restore_start(hybrid, incoming, budget, *, positive=False):
    """Hold both graph inequalities hard during the shared SLSQP Phase I."""
    scaled = UnitBoxEvaluator(hybrid)
    if exact_feasible(hybrid, incoming, positive=positive):
        return incoming, []
    attempts = []
    for margin in (.05, 0.):
        outcome, cost = measured(lambda: restore(scaled, incoming, margin=margin,
                                                 max_evaluations=budget, hard_tail=2), hybrid)
        feasible_designs = [outcome.accepted_design]
        feasible_designs.extend(point.design for point in outcome.trial_history
                                if point.engineering_feasible and point.margin_feasible)
        selected = next((np.asarray(design) for design in feasible_designs
                         if design is not None and exact_feasible(hybrid, design, positive=positive)), None)
        attempts.append({"result": asdict(outcome), "cost": cost, "selected_unit_design": selected})
        if selected is not None:
            return selected, attempts
    return None, attempts


def shared_seed(output, evaluator, model, budget):
    """Use one positive baffle-only start and restore through the common constraints."""
    records = baffle_sweep(evaluator, model)
    preserve(output, "baffle_only_sweep", records)
    eligible = [row for row in records if row["height_mm"] > MIN_SEED_HEIGHT_MM
                and row.get("reconstructed_connected") and all(row.get("strict_graph_paths", (False,)))
                and row.get("fitted_sharp_baffle", {}).get("state_passed")
                and row["reconstructed_area_mm2"] <= MAX_AREA_MM2
                and max(row["hybrid_immersed"]["inequalities"][2:]) <= FEASIBILITY_TOLERANCE]
    fitted_order = sorted(eligible, key=lambda row: abs(row["fitted_sharp_baffle"]["side_share"] - .5))
    immersed_order = sorted(eligible, key=lambda row: abs(row["hybrid_immersed"]["side_fraction"] - .5))
    feasible_order = [row for row in immersed_order
                      if max(row["hybrid_immersed"]["inequalities"]) <= FEASIBILITY_TOLERANCE]
    selected_heights, ranked_starts = set(), []
    for group in (feasible_order[:1], immersed_order[:1], fitted_order[:1], immersed_order, fitted_order):
        for row in group:
            if row["height_mm"] not in selected_heights:
                selected_heights.add(row["height_mm"])
                ranked_starts.append(row)
    records = []
    selected = None
    for row in ranked_starts[:3]:
        incoming = unit_design(row["physical_design"])
        selected, attempts = restore_start(evaluator, incoming, budget, positive=True)
        records.append({"initial_baffle_height_mm": row["height_mm"], "restoration": attempts,
                        "restored_unit_design": selected})
        if selected is not None:
            break
    preserve(output, "shared_restoration", records)
    return {"shared_unit_design": selected, "shared_physical_design": None if selected is None else physical_design(selected),
            "successful": selected is not None, "attempted_initial_heights_mm": [row["initial_baffle_height_mm"] for row in records]}


def audit_candidate(hybrid, scaled_design, output):
    scaled = UnitBoxEvaluator(hybrid)
    response, cost = measured(lambda: scaled.evaluate(scaled_design), hybrid)
    design = physical_design(scaled_design)
    height_mm = HEIGHT_MM * design[-1]
    density = hybrid.density(design)
    certificate = hybrid.certificate(design)
    dimensional_power, fluxes = _power_and_flux(hybrid.domain, hybrid.state.values)
    fraction = float(fluxes[2] / fluxes[0])
    objective = dimensional_power / hybrid.power_reference + hybrid.height_cost * design[-1]
    np.testing.assert_allclose([dimensional_power, fraction, objective],
                               [response.power_w_per_m, response.side_fraction, response.objective], rtol=1e-10, atol=1e-12)
    inequalities = response.constraints.copy()
    inequalities[:2] = [(0.499 - fraction) / .001, (fraction - 0.501) / .001]
    augmented_kkt = audit_kkt(scaled_design, response.gradient, inequalities, response.jacobian)
    hydraulic_kkt = audit_kkt(scaled_design, response.gradient, inequalities[:3], response.jacobian[:3])
    geometry = reconstruct_hybrid(hybrid.pipeline.grid, density.projected, height_mm)
    contour = write_reconstruction(geometry, output / "geometry.json")
    feasible = (max(inequalities) <= FEASIBILITY_TOLERANCE and all(certificate.hard_paths)
                and geometry.port_component is not None and geometry.area_mm2 <= MAX_AREA_MM2)
    return dict(response_summary(response), unit_design=scaled_design, physical_design=design, baffle_height_mm=height_mm,
                baffle_height_ratio=design[-1], filtered_density=density.filtered, projected_density=density.projected,
                primary_hydraulic_kkt=asdict(hydraulic_kkt), augmented_kkt=asdict(augmented_kkt),
                independently_integrated_side_fraction=fraction, independently_integrated_power_W_per_m=dimensional_power,
                independently_integrated_objective=objective,
                exact_and_immersed_feasible=feasible, graph_conductances=certificate.conductance,
                strict_graph_paths=certificate.hard_paths, geometry=contour, cost=cost), geometry


def run_model(output, context, model_name, initial, stages, mma_budget, restore_budget):
    domain, model, reference, initial_state, graph = context
    selected, state = initial, initial_state
    summaries = []
    for stage_number, stage in enumerate(stages):
        folder = output / f"stage_{stage_number:02d}"
        folder.mkdir()
        evaluator = HybridEvaluator(domain, model, stage, reference, initial_state=state,
                                    height_cost=HEIGHT_COSTS[model_name], graph=graph)
        start, attempts = restore_start(evaluator, selected, restore_budget)
        preserve(folder, "restoration", {"stage": asdict(stage), "initial_unit_design": selected,
                                         "restored_unit_design": start, "attempts": attempts})
        if start is None:
            summaries.append({"stage": asdict(stage), "reason": "no feasible stage handoff"})
            break
        wrapper = UnitBoxEvaluator(evaluator)
        outcome, cost = measured(lambda: optimize(wrapper, start, algorithm="mma", max_evaluations=mma_budget), evaluator)
        preserve(folder, "mma", {"stage": asdict(stage), "result": asdict(outcome), "cost": cost})
        candidates = []
        for name, design in (("terminal", outcome.terminal_design),
                             ("best_trial", outcome.best_feasible_design), ("restored_start", start)):
            if design is None or not np.isfinite(design).all() or np.any((design < 0.) | (design > 1.)):
                continue
            if any(np.array_equal(design, row[1]) for row in candidates):
                continue
            audit_dir = folder / name
            audit_dir.mkdir()
            independent = HybridEvaluator(domain, model, stage, reference, initial_state=state,
                                          height_cost=HEIGHT_COSTS[model_name], graph=graph)
            try:
                audit, geometry = audit_candidate(independent, design, audit_dir)
            except (RuntimeError, ValueError, np.linalg.LinAlgError) as error:
                preserve(audit_dir, "checkpoint", {"unit_design": design, "error": f"{type(error).__name__}: {error}"})
                continue
            preserve(audit_dir, "checkpoint", audit)
            candidates.append((name, design, audit, geometry, independent.state))
        feasible = [row for row in candidates if row[2]["exact_and_immersed_feasible"]]
        if not feasible:
            summaries.append({"stage": asdict(stage), "reason": "no independently confirmed feasible endpoint"})
            break
        stationary = [row for row in feasible if row[2]["augmented_kkt"]["verified"]]
        name, selected, checkpoint, geometry, state = min(stationary or feasible, key=lambda row: row[2]["objective"])
        fitted = fitted_screen(geometry)
        report = {"stage": asdict(stage), "selected_source": name, "selected": checkpoint,
                  "physical_fitted_screen": fitted,
                  "candidate_directories": [row[0] for row in candidates]}
        preserve(folder, "selected", report)
        summaries.append(report)

    finalists = [row for row in summaries if "selected" in row]
    last = finalists[-1] if finalists else None
    if last is not None:
        archive = output / "final_export"
        archive.mkdir()
        checkpoint = last["selected"]
        export = export_endpoint(archive, grid=model.grid, projected=checkpoint["projected_density"],
                                 baffle_height_mm=checkpoint["baffle_height_mm"], domain=domain, state=state,
                                 fitted_spacings_mm=FITTED_SPACINGS_MM)
        preserve(archive, "export_report", export)
    return {"model": model_name, "height_cost": HEIGHT_COSTS[model_name], "stage_results": summaries,
            "all_stages_selected": len(finalists) == len(stages), "last_selected": last}


def run(output, route, max_stages, mode, shared_seed_report, mma_budget, restore_budget):
    grid = DesignGrid(.5)
    graph = ConnectivityGraph(grid, baffle=RectangularBaffleClearance(
        x_min_mm=BAFFLE_X_MM[0], x_max_mm=BAFFLE_X_MM[1], max_height_mm=HEIGHT_MM * HEIGHT_LIMIT))
    domain = build_domain(grid, .5)
    model = HybridDensityModel(grid)
    reference, state = hybrid_reference(domain, model, BETA_ZERO_STAGE.alpha_max)
    raw_reference_power, raw_reference_flux = _power_and_flux(domain, state.values)
    evaluator = HybridEvaluator(domain, model, BETA_ZERO_STAGE, reference, initial_state=state, graph=graph)
    zero_height = evaluator.evaluate(model.initial)
    stages = [Stage(.01, beta, alpha) for beta, alpha in ROUTES[route]][:max_stages]
    setup = {"slurm_job_id": os.environ["SLURM_JOB_ID"], "design_spacing_mm": grid.spacing_mm,
             "mesh_triangles": len(domain.fem.triangles), "design_variables": model.size,
             "stages": [asdict(stage) for stage in stages], "physical_zero_height_design": model.initial,
             "normalized_power_reference_W_per_m": reference,
             "legacy_sharp_zero_height_reference_W_per_m": raw_reference_power,
             "legacy_sharp_zero_height_side_fraction": float(raw_reference_flux[2] / raw_reference_flux[0]),
             "hybrid_sharp_zero_height_response": response_summary(zero_height),
             "height_coordinate_units": "s=h/4 mm; unit-box coordinate z=s/0.9375",
             "height_costs_by_model": HEIGHT_COSTS,
             "mma_budget_per_stage": mma_budget, "restoration_budget_per_attempt": restore_budget,
             "graph_threshold": graph.threshold, "graph_cutoff": graph.cutoff,
             "graph_floor_budget": graph.floor_budget, "graph_minimum_conductance": graph.minimum_conductance}
    preserve(output, "setup", setup)
    if shared_seed_report is None:
        seed = shared_seed(output, evaluator, model, restore_budget)
        preserve(output, "seed", seed)
    else:
        stored = json.loads(shared_seed_report.read_text(encoding="utf-8"))
        if not stored["completed"]:
            raise ValueError("The archived common seed did not run successfully")
        seed = stored["results"]
    result = {"diagnostic_only": True, "common_seed": seed, "model_results": {}}
    if mode == "seed" or not seed["successful"]:
        return result
    common = np.asarray(seed["shared_unit_design"], dtype=float)
    if common.shape != (model.size,) or not np.isfinite(common).all() or np.any((common < 0.) | (common > 1.)):
        raise ValueError("Archived common seed has invalid unit-box coordinates")
    if not exact_feasible(evaluator, common, positive=True):
        raise ValueError("Archived common seed does not reproduce a positive feasible hybrid design")
    context = domain, model, reference, evaluator.state, graph
    for model_name in (tuple(HEIGHT_COSTS) if mode == "all" else (mode,)):
        lane = output / model_name
        lane.mkdir()
        completed = write_report(lane, "summary", lambda: run_model(
            lane, context, model_name, common.copy(), stages, mma_budget, restore_budget))
        result["model_results"][model_name] = {"completed": completed, "report": f"{model_name}/summary.json"}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="New per-job evidence directory")
    parser.add_argument("--mode", choices=("seed", "model-a", "model-b", "all"), default="all")
    parser.add_argument("--route", choices=ROUTES, default="coupled")
    parser.add_argument("--max-stages", type=int)
    parser.add_argument("--shared-seed", type=Path, help="Existing completed seed.json from the common diagnostic")
    parser.add_argument("--mma-evaluations", type=int, default=150)
    parser.add_argument("--restoration-evaluations", type=int, default=250)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical diagnostic as a Sol Slurm job")
    if min(args.mma_evaluations, args.restoration_evaluations) < 1 or (args.max_stages is not None and args.max_stages < 1):
        parser.error("Stage and evaluation budgets must be positive")
    if args.mode == "seed" and args.shared_seed is not None:
        parser.error("A seed-only run must compute its own seed")
    if args.mode in HEIGHT_COSTS and args.shared_seed is None:
        parser.error("A single-model run requires the common seed.json")
    args.output.mkdir(parents=True, exist_ok=False)
    preserve(args.output, "environment", environment())
    if not write_report(args.output, "summary", lambda: run(
            args.output, args.route, args.max_stages, args.mode, args.shared_seed,
            args.mma_evaluations, args.restoration_evaluations)):
        raise SystemExit(1)
    report = json.loads((args.output / "summary.json").read_text(encoding="utf-8"))["results"]
    if any(not lane["completed"] for lane in report["model_results"].values()):
        raise SystemExit(1)
    if not report["common_seed"]["successful"]:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
