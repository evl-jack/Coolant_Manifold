"""Independently verify the zero-inertia fixed-pipe state and responses."""

import argparse
import json
from pathlib import Path

import numpy as np

from ..evaluator import Stage
from ..stokes_fixed_pipe_model import StokesFixedPipeEvaluator, build_stokes_fixed_pipe
from .fixed_pipe_adjoint_audit import (
    STEP_SIZES, directional_check, filter_transfer_checks, independent_filter, state_jacobian_checks,
)


STOKES_STEPS = tuple(sorted((*STEP_SIZES, 5e-4, 2e-4), reverse=True))


def audit(pilot, selected=None):
    archive_case = pilot / "stage_03_beta_8" / "original"
    reference = json.loads((archive_case / "metrics.json").read_text(encoding="utf-8"))["results"]
    domain, pipeline = build_stokes_fixed_pipe()
    evaluator = StokesFixedPipeEvaluator(domain, pipeline, Stage(.10, 8., 1e5),
                                        reference["power_w_per_m"], reference["side_fraction"])
    rng = np.random.default_rng(473208)
    original = pipeline.initial
    with np.load(archive_case / "fields.npz", allow_pickle=False) as archive:
        np.testing.assert_array_equal(original, archive["design"])
    gray = .15 + .70 * original + rng.uniform(-.02, .02, pipeline.size)
    cases = {"original": original, "interior_gray": gray}
    if selected is not None:
        with np.load(selected, allow_pickle=False) as snapshot:
            cases["selected"] = snapshot["design"]

    matrix, fixed_rhs = independent_filter(pipeline)
    filter_checks = filter_transfer_checks(evaluator, matrix, fixed_rhs, gray, rng)
    report = {"fem_advection_coefficient": domain.fem.inertia, "power_reference_W_per_m": evaluator.power_reference,
              "side_reference": evaluator.side_reference, "filter": filter_checks, "cases": {}}
    for label, design in cases.items():
        value = evaluator.evaluate(design)
        state = evaluator.state
        resistance = evaluator.transfer.forward(evaluator.density(design).resistance)
        state_check = state_jacobian_checks(evaluator, state, resistance, rng)
        no_convection = all(all(norm == 0. for term, norm in row["term_inf_norms"].items()
                                if term.startswith("convection"))
                            for row in state_check["directions"].values())
        tangent = directional_check(evaluator, matrix, design, value, state,
                                    rng.uniform(-1, 1, pipeline.size), STOKES_STEPS)
        report["cases"][label] = {"state_jacobian": state_check, "zero_convection_actions": no_convection,
                                  "power_W_per_m": value.power_w_per_m, "side_fraction": value.side_fraction,
                                  "residual": value.nonlinear_residual, "mass_imbalance": value.mass_imbalance,
                                  "linear_updates": len(state.residual_history) - 1, "directional_check": tangent,
                                  "passed": bool(state_check["passed"] and no_convection and tangent["passed"])}
        print(f"{label}: linear state {state_check['passed']}, derivative {tangent['derivative_passed']}, "
              f"Taylor {tangent['taylor_passed']}", flush=True)
    report["passed"] = bool(filter_checks["passed"] and all(row["passed"] for row in report["cases"].values()))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pilot", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--selected", type=Path)
    args = parser.parse_args()
    report = audit(args.pilot, args.selected)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    raise SystemExit(not report["passed"])


if __name__ == "__main__":
    main()
