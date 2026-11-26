"""Search for a connected full-density feasibility seed, with no production changes."""

import argparse
from dataclasses import asdict, replace
import os
from pathlib import Path

import numpy as np

from .connectivity import ConnectivityDiagnosticEvaluator
from .density import DesignGrid, DensityPipeline
from .domain import build_domain
from .evaluator import Evaluator, Stage, baseline_reference
from .expanded_envelope import ExpandedDesignGrid, build_expanded_domain
from .full_density_continuation import evaluate_candidate
from .gate2 import write_report
from .optimize import optimize


class SplitRestoration:
    """Minimize distance to the midpoint of the original split band."""

    def __init__(self, evaluator):
        self.evaluator = evaluator

    def evaluate(self, design):
        response = self.evaluator.evaluate(design)
        difference = response.side_fraction - 0.5
        fraction_gradient = 0.001 * response.jacobian[1]
        return replace(response, objective=difference**2, gradient=2 * difference * fraction_gradient,
                       constraints=response.constraints[2:], jacobian=response.jacobian[2:])


def run(output, envelope):
    if envelope == "expanded":
        grid = ExpandedDesignGrid(0.5)
        domain = build_expanded_domain(grid, 0.5)
    else:
        grid = DesignGrid(0.5)
        domain = build_domain(grid, 0.5)
    pipeline = DensityPipeline(grid)
    alpha = 1e6
    reference, seed = baseline_reference(domain, pipeline, alpha)
    stage = Stage(0.01, 0., alpha)
    hydraulic = Evaluator(domain, pipeline, stage, reference, initial_state=seed)
    augmented = ConnectivityDiagnosticEvaluator(hydraulic)
    restoration = SplitRestoration(augmented)
    initial = augmented.evaluate(grid.base_raw)
    outcome = optimize(restoration, grid.base_raw, algorithm="slsqp", max_evaluations=250)
    if not write_report(output, "split_search", lambda: asdict(outcome)):
        raise RuntimeError("Could not archive split-search history")
    selected = outcome.best_feasible_design
    result = dict(envelope=envelope, stage=asdict(stage), initial_inequalities=initial.constraints,
                  initial_side_fraction=initial.side_fraction, diagnostic_objective="(side_fraction - 0.5)^2",
                  physical_feasible_candidate=False, phase_two_executed=False,
                  slurm_job_id=os.environ["SLURM_JOB_ID"])
    if selected is not None:
        candidate = output / "best_area_and_graph_trial"
        candidate.mkdir()
        checked, _ = evaluate_candidate(domain, pipeline, stage, reference, seed, selected, candidate,
                                        evaluate_fitted=True, connectivity_diagnostic=True)
        if not write_report(candidate, "checkpoint", lambda: checked):
            raise RuntimeError("Could not archive independently replayed split candidate")
        result.update(best_immersed_side_fraction=checked["side_fraction"],
                      original_engineering_feasible=checked["engineering_feasible"],
                      connected_ports=checked.get("reconstruction", {}).get("connected_ports"),
                      reconstructed_area_mm2=checked.get("reconstruction", {}).get("physical_fluid_area_mm2"),
                      fitted_screen=checked.get("fitted_screen"),
                      physical_feasible_candidate=bool(
                          checked["engineering_feasible"] and checked.get("connectivity_diagnostic", {}).get(
                              "connectivity_only_satisfied") and checked.get("fitted_screen", {}).get(
                                  "fitted_diagnostic_passed")))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("envelope", choices=("original", "expanded"))
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical diagnostic on Sol")
    args.output.mkdir(parents=True)
    if not write_report(args.output, "summary", lambda: run(args.output, args.envelope)):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
