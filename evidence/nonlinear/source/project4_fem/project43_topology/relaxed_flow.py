"""Matched full-density power minimization on the original T-junction."""

import argparse
from dataclasses import asdict, replace
import json
import os
from pathlib import Path

import numpy as np

from ..fem import State
from .connectivity import ConnectivityDiagnosticEvaluator, ConnectivityGraph
from .density import DensityPipeline, DesignGrid
from .domain import build_domain
from .evaluator import Evaluator, Stage, baseline_reference
from .full_density_continuation import ROUTES, fitted_screen, restore_stage
from .gate2 import environment, write_report
from .gate2_restoration import FEASIBILITY_TOLERANCE, measured
from .kkt_audit import audit_kkt
from .optimize import optimize
from .reconstruction import Reconstruction, reconstruct, write_reconstruction
from .verify_flow import _power_and_flux


BANDS = {"model-a": (.20, .80), "model-b": (.10, .90)}
SPLIT_SCALE = .001
MMA_BUDGET = 150
RESTORE_BUDGET = 250


class BandedEvaluator(Evaluator):
    """Keep the existing derivatives and change only the fraction bounds."""

    def __init__(self, domain, pipeline, stage, reference, *, band, initial_state=None):
        super().__init__(domain, pipeline, stage, reference, initial_state=initial_state)
        self.band = band

    def evaluate(self, design):
        value = super().evaluate(design)
        lower, upper = self.band
        fraction = value.side_fraction
        constraints = np.r_[(lower - fraction) / SPLIT_SCALE, (fraction - upper) / SPLIT_SCALE,
                            value.constraints[2]]
        return replace(value, constraints=constraints)


def save(output, name, data):
    if not write_report(output, name, lambda: data):
        raise RuntimeError(f"Could not archive {name} in {output}")


def setup():
    grid = DesignGrid(.5)
    return build_domain(grid, .5), DensityPipeline(grid), ConnectivityGraph(grid)


def prepare(output):
    from shapely.geometry import box

    domain, pipeline, _ = setup()
    power, seed = baseline_reference(domain, pipeline, 1e5)
    np.savez(output / "baseline.npz", power=power, state=seed.values)
    polygon = box(0., 0., 20., 4.).union(box(9.5, 4., 10.5, 7.))
    geometry = Reconstruction(polygon, (polygon,), (0, 0, 0), ())
    save(output, "sharp_fitted_baseline",
         fitted_screen(geometry, (0., 1.), (.25, .125, .0625)))
    return {"normalization_W_per_m": power, "design_variables": pipeline.grid.size,
            "triangles": len(domain.fem.triangles), "area_cap_mm2": 83.}


def independent_replay(hydraulics, graph, design, output, reference_fitted_power):
    augmented = ConnectivityDiagnosticEvaluator(hydraulics, graph)
    response = augmented.evaluate(design)
    density = hydraulics.density(design)
    power, flux = _power_and_flux(hydraulics.domain, hydraulics.state.values)
    fraction = float(flux[2] / flux[0])
    area = 41. + hydraulics.pipeline.grid.cell_area_mm2 * density.projected.sum()
    low, high = hydraulics.band
    constraints = np.r_[(low - fraction) / SPLIT_SCALE, (fraction - high) / SPLIT_SCALE,
                        area - 83., response.constraints[3:]]
    objective = power / hydraulics.power_reference
    np.testing.assert_allclose(objective, response.objective, rtol=1e-10, atol=1e-12)
    np.testing.assert_allclose(constraints, response.constraints, rtol=1e-10, atol=1e-9)
    certificate, _ = augmented.certificate(design)
    kkt = audit_kkt(design, response.gradient, constraints, response.jacobian)
    geometry = reconstruct(hydraulics.pipeline.grid, density.projected)
    contour = write_reconstruction(geometry, output / "geometry.json")
    admissible = kkt.primal_violation <= FEASIBILITY_TOLERANCE and all(certificate.hard_paths)
    screen = fitted_screen(geometry, hydraulics.band) if admissible else {}
    finest = screen.get("levels", [])
    fitted_power = finest[-1]["mechanical_power_W_per_m"] if screen.get("completed") else None
    return {"design": design, "objective": objective, "immersed_power_W_per_m": power,
            "immersed_side_share": fraction, "immersed_area_mm2": area,
            "inequalities": constraints, "strict_graph_paths": certificate.hard_paths,
            "immersed_feasible": admissible, "kkt": asdict(kkt), "geometry": contour,
            "fitted": screen, "fitted_passed": bool(screen.get("fitted_diagnostic_passed")),
            "fitted_power_W_per_m": fitted_power,
            "fitted_power_ratio_vs_sharp": None if fitted_power is None else fitted_power / reference_fitted_power}


def candidate_designs(outcome, start):
    ordered = sorted(outcome.feasible_history, key=lambda point: point.objective)
    # Keep one economical point and several with margin against reference-model bias.
    margin = sorted(ordered, key=lambda point: (max(point.inequalities[:2]), point.objective))
    options = [("terminal", outcome.terminal_design)]
    if ordered:
        options.append(("best_feasible", ordered[0].design))
    if margin:
        options.extend([("split_margin", margin[0].design), ("intermediate", margin[len(margin) // 2].design)])
    options.append(("restored_start", start))
    unique = []
    for label, design in options:
        if not any(np.array_equal(design, other) for _, other in unique):
            unique.append((label, design))
    return unique


def run(output, model, common, *, inertia=50.):
    from Project43_Accelerated.hybrid import accelerated_domain

    domain, pipeline, graph = setup()
    wrapped = accelerated_domain(domain, device="cpu")
    domain.fem.inertia = wrapped.fem.inertia = inertia
    band = BANDS[model]
    with np.load(common / "baseline.npz") as data:
        reference = float(data["power"])
        values = data["state"].copy()
    free = np.setdiff1d(np.arange(domain.fem.size), domain.fixed)
    seed = State(values, [], None, free)
    fitted_baseline = json.loads((common / "sharp_fitted_baseline.json").read_text())["results"]
    if not fitted_baseline["fitted_diagnostic_passed"]:
        raise RuntimeError("The independent sharp fitted baseline failed its mesh or state screen")
    reference_fitted_power = fitted_baseline["levels"][-1]["mechanical_power_W_per_m"]
    incoming = pipeline.grid.base_raw.copy()
    stages = [Stage(.01, beta, alpha) for beta, alpha in ROUTES["coupled"]]
    save(output, "setup", {"model": model, "flow_band": band, "inertia": inertia,
         "design_cells": pipeline.grid.size, "area_cap_mm2": 83., "mma_budget": MMA_BUDGET,
         "restoration_budget": RESTORE_BUDGET, "stages": list(map(asdict, stages)),
         "normalization_W_per_m": reference, "fitted_sharp_power_W_per_m": reference_fitted_power})
    best = None
    completed = []
    for index, stage in enumerate(stages):
        stage_dir = output / f"stage_{index:02d}"
        stage_dir.mkdir()
        hydraulic = BandedEvaluator(wrapped, pipeline, stage, reference, band=band, initial_state=seed)
        evaluator = ConnectivityDiagnosticEvaluator(hydraulic, graph)
        start, initial, restoration = restore_stage(evaluator, hydraulic, incoming, RESTORE_BUDGET, True)
        save(stage_dir, "restoration", {"accepted_design": start, "initial": initial, "attempts": restoration})
        if start is None:
            break
        outcome, work = measured(lambda: optimize(evaluator, start, algorithm="mma",
                                                 max_evaluations=MMA_BUDGET), hydraulic)
        save(stage_dir, "mma", {"result": asdict(outcome), "cost": work})
        audited = []
        for label, design in candidate_designs(outcome, start):
            candidate_dir = stage_dir / label
            candidate_dir.mkdir()
            replay = BandedEvaluator(domain, pipeline, stage, reference, band=band, initial_state=seed)
            try:
                report = independent_replay(replay, graph, design, candidate_dir, reference_fitted_power)
            except (RuntimeError, ValueError, AssertionError, np.linalg.LinAlgError) as error:
                save(candidate_dir, "checkpoint", {"error": f"{type(error).__name__}: {error}"})
                continue
            save(candidate_dir, "checkpoint", report)
            relative = str(candidate_dir.relative_to(output))
            if report["fitted_passed"]:
                qualified = {"checkpoint": relative, "fitted_power_W_per_m": report["fitted_power_W_per_m"],
                             "fitted_power_ratio_vs_sharp": report["fitted_power_ratio_vs_sharp"],
                             "fitted_side_share": report["fitted"]["levels"][-1]["side_share"],
                             "fitted_area_mm2": report["geometry"]["physical_fluid_area_mm2"],
                             "optimizer_sample": not np.array_equal(design, start),
                             "kkt_passed": report["kkt"]["verified"]}
                if best is None or (not qualified["optimizer_sample"], qualified["fitted_power_W_per_m"]) < (
                        not best["optimizer_sample"], best["fitted_power_W_per_m"]):
                    best = qualified
                save(candidate_dir, "physical_acceptance", qualified)
            if (report["immersed_feasible"] and report["geometry"]["connected_ports"]
                    and report["geometry"]["physical_fluid_area_mm2"] <= 83.):
                audited.append((report, replay.state, relative))
        if not audited:
            break
        report, seed, path = min(audited, key=lambda row: (not row[0]["fitted_passed"], row[0]["objective"]))
        incoming = report["design"]
        completed.append(path)
        save(stage_dir, "selected", {"checkpoint": path, "best_physical_so_far": best})
        print(f"{model} stage {index} complete, best physical {best}", flush=True)
    return {"completed_stage_selections": completed, "scheduled_stages": len(stages),
            "best_physical": best, "model": model}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("prepare", *BANDS))
    parser.add_argument("output", type=Path)
    parser.add_argument("--common", type=Path)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run on Sol")
    if args.mode != "prepare" and args.common is None:
        parser.error("Supply the completed common baseline directory")
    args.output.mkdir(parents=True)
    write_report(args.output, "environment", environment)
    if not write_report(args.output, "summary", lambda: (
            prepare(args.output) if args.mode == "prepare" else run(args.output, args.mode, args.common))):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
