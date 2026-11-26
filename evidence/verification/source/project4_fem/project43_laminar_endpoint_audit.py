"""Audit the archived Sol 2350992 CFD endpoint in its original SLSQP units."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import nnls


SOURCE_HASH = "e409c6c4fbc39de2721d6a73241dfac965c4a536b61e2d41f94c5b4b3f57b40e"
RESULT_HASH = "d01f778c24a668f1e76a9a656d052c29720dc1b70fb008b8c0e65cd82ef82578"
ACTIVE_TOLERANCE = 1e-6
GATES = dict(primal_violation=1e-7, dual_violation=1e-7, stationarity_inf=1e-5,
             complementarity_inf=1e-6)


def verified_bytes(path, expected):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f"Unexpected archived checksum for {path}")
    return data


def archived_evaluator(source, design):
    # Use the exact archived solver, discretization and response derivative.
    # Return its evaluator before the old driver starts SLSQP.
    data = verified_bytes(source, SOURCE_HASH).decode()
    marker = "    start=np.array([throat,width]); r0,j0=evaluate(start)"
    setup, matched, _ = data.partition(marker)
    if not matched:
        raise ValueError("Cannot isolate the archived response evaluator")
    namespace = {"__name__": "project4_fem.archived_lam_pair", "__package__": "project4_fem"}
    exec(compile(setup + "    return evaluate\n", str(source), "exec"), namespace)
    return namespace["run"]({"best_feasible": dict(open_height_mm=design[0], shoulder_width_mm=design[1])})


def audit(archive):
    historical = json.loads(verified_bytes(archive / "pair.json", RESULT_HASH))
    x = np.array([historical["h"], historical["s"]])
    evaluate = archived_evaluator(archive / "project4_fem/lam_pair.py", x)
    value, derivative = evaluate(x)
    np.testing.assert_allclose(value, historical["density"], rtol=2e-8, atol=1e-12)
    steps = []
    for step in (1e-3, 5e-4):
        observed = np.column_stack([(evaluate(x + delta)[0] - evaluate(x - delta)[0])/(2*step)
                                    for delta in step*np.eye(2)])
        relative = np.abs(observed-derivative)/np.maximum(np.abs(derivative), 1e-12)
        steps.append(dict(step_mm=step, gradient=observed.tolist(), relative_error=relative.tolist(),
                          passed=bool(np.all(relative <= 1e-3))))

    reference = historical["start"][0]
    gradient = derivative[0]/reference
    h, s = x
    q = value[1]
    names = ("split_lower", "split_upper", "area_cap", "throat_lower", "shoulder_lower",
             "throat_upper", "shoulder_upper")
    constraints = np.array([100*(q-.499), 100*(.501-q), 4-h-4*s, h-.15, s, 4-h, 1-s])
    jacobian = np.vstack((100*derivative[1], -100*derivative[1], [-1., -4.], np.eye(2), -np.eye(2)))
    active = constraints <= ACTIVE_TOLERANCE
    active_jacobian = jacobian[active]
    multipliers = np.zeros(len(names))
    multipliers[active], _ = nnls(active_jacobian.T, gradient)
    stationarity = gradient - multipliers @ jacobian
    numerical = dict(primal_violation=max(0., -float(constraints.min())),
                     dual_violation=max(0., -float(multipliers.min())),
                     stationarity_inf=float(np.linalg.norm(stationarity, np.inf)),
                     complementarity_inf=float(np.linalg.norm(multipliers*constraints, np.inf)))
    rank = int(np.linalg.matrix_rank(active_jacobian))
    constraints_report = {name: dict(slack=float(c), gradient=g.tolist(), active=bool(a), multiplier=float(lam))
                          for name, c, g, a, lam in zip(names, constraints, jacobian, active, multipliers)}
    return dict(archived_sol_job=2350992, design_mm=x.tolist(), reference_power_W_m=reference,
                power_W_m=float(value[0]), side_share_fraction=float(q), nominal_area_mm2=float(79+h+4*s),
                objective_gradient_per_mm=gradient.tolist(), raw_response_gradients=derivative.tolist(),
                endpoint_replay_difference=(value-historical["density"]).tolist(),
                endpoint_gradient_checks=steps, historical_termination=historical["message"],
                historical_iterations=historical["iterations"], active_tolerance=ACTIVE_TOLERANCE,
                constraints=constraints_report, active_gradient_rank=rank,
                independent_active_gradients=rank == int(active.sum()),
                active_gradient_condition_number=float(np.linalg.cond(active_jacobian)),
                stationarity_vector=stationarity.tolist(), residuals=numerical, gates=GATES,
                passed=all(step["passed"] for step in steps)
                and rank == int(active.sum()) and all(numerical[k] <= gate for k, gate in GATES.items()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="Frozen lam-pair-01 directory from Sol")
    args = parser.parse_args()
    result = audit(args.archive)
    print(json.dumps(result, indent=2, allow_nan=False))
    raise SystemExit(not result["passed"])
