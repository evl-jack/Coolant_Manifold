"""Independent first-order KKT audit for the two-control inequalities."""

import numpy as np
from scipy.optimize import nnls

from .project43_optimizer import BOUNDS


def kkt_report(evaluator, x, result=None, active_tolerance=1e-6):
    value = evaluator.evaluate(x)
    x = np.asarray(x, dtype=float)
    lower = x - np.array([bound[0] for bound in BOUNDS])
    upper = np.array([bound[1] for bound in BOUNDS]) - x
    inequalities = np.r_[value.constraints, lower, upper]
    gradients = np.vstack((value.jacobian, np.eye(2), -np.eye(2)))
    names = ["side_share_lower", "side_share_upper", "nominal_area",
             "x1_lower", "x2_lower", "x1_upper", "x2_upper"]
    active = inequalities <= active_tolerance
    multipliers = np.zeros(len(names))
    if np.any(active):
        multipliers[active], _ = nnls(gradients[active].T, value.gradient)
    stationarity = value.gradient - multipliers @ gradients
    primal = max(0., -float(np.min(inequalities)))
    dual = max(0., -float(np.min(multipliers)))
    report = {
        "x": x.tolist(), "objective": value.objective,
        "side_share": value.side_share, "nominal_area_mm2": float(79 + 4*x.sum()),
        "active_tolerance": active_tolerance,
        "active_constraints": [name for name, flag in zip(names, active) if flag],
        "multipliers": dict(zip(names, multipliers.tolist())),
        "primal_violation": primal, "dual_violation": dual,
        "stationarity_inf": float(np.max(np.abs(stationarity))),
        "complementarity_inf": float(np.max(np.abs(multipliers * inequalities))),
        "slsqp_success": None if result is None else bool(result.success),
        "slsqp_message": None if result is None else str(result.message),
        "unique_design_evaluations": evaluator.unique_evaluations,
    }
    report["kkt_verified"] = (primal <= 1e-7 and dual <= 1e-7 and
                              report["stationarity_inf"] <= 1e-5 and
                              report["complementarity_inf"] <= 1e-6)
    return report
