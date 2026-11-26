"""Independent first-order audit for unit-box designs with g <= 0."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import linprog


@dataclass(frozen=True)
class KKTReport:
    verified: bool
    multipliers: np.ndarray
    multiplier_source: str
    lp_status: int | None
    lp_message: str | None
    inequality_violation: float
    box_violation: float
    primal_violation: float
    dual_violation: float
    stationarity_inf: float
    complementarity_inf: float
    projected_residual: np.ndarray


def audit_kkt(
    design,
    objective_gradient,
    inequalities,
    constraint_jacobian,
    *,
    primal_tolerance=1e-6,
    dual_tolerance=1e-8,
    stationarity_tolerance=1e-4,
    complementarity_tolerance=1e-5,
):
    """Recover a multiplier witness without using the optimizer's multipliers.

    With d = df + Jg.T @ lambda, the signed projected residual at each
    coordinate is clip(d_i, design_i - 1, design_i). For an in-box design,
    the stationarity bound is therefore equivalent to an upper bound on
    d_i when design_i > tolerance and a lower bound when 1 - design_i >
    tolerance. Construct those bounds at the actual design coordinates.
    """
    design = np.asarray(design, dtype=float)
    gradient = np.asarray(objective_gradient, dtype=float)
    inequalities = np.asarray(inequalities, dtype=float)
    jacobian = np.asarray(constraint_jacobian, dtype=float)
    tolerance = stationarity_tolerance
    tolerances = [primal_tolerance, dual_tolerance, tolerance, complementarity_tolerance]
    if not np.all(np.isfinite(tolerances)) or min(tolerances) < 0:
        raise ValueError("KKT tolerances must be finite and nonnegative")
    if (
        design.ndim != 1
        or not design.size
        or gradient.shape != design.shape
        or inequalities.ndim != 1
        or jacobian.shape != (inequalities.size, design.size)
        or not all(np.all(np.isfinite(a)) for a in (design, gradient, inequalities, jacobian))
    ):
        raise ValueError("KKT inputs must have consistent shapes and finite values")

    inequality_violation = float(np.max(inequalities, initial=0.0))
    box_violation = float(max(0.0, -design.min(), design.max() - 1.0))
    primal = max(inequality_violation, box_violation)

    multipliers = np.zeros(inequalities.size)
    source = "zero multipliers"
    lp_status = None
    lp_message = None
    projected = design - np.clip(design - gradient, 0.0, 1.0)
    if inequalities.size and np.max(np.abs(projected)) > tolerance:
        upper = design > tolerance
        lower = 1.0 - design > tolerance
        coefficients = np.vstack((jacobian[:, upper].T, -jacobian[:, lower].T))
        lhs = np.column_stack((coefficients, -np.ones(len(coefficients))))
        rhs = np.concatenate((-gradient[upper], gradient[lower]))

        # Translate |lambda_j * g_j| <= complementarity_tolerance into
        # multipliers' bounds without treating near-zero constraints as zero.
        absolute_g = np.abs(inequalities)
        caps = np.full(inequalities.size, np.inf)
        with np.errstate(over="ignore"):
            np.divide(complementarity_tolerance, absolute_g, out=caps, where=absolute_g > 0)
        bounds = [(0.0, float(cap) if np.isfinite(cap) else None) for cap in caps]
        # Seek an interior witness when possible. Minimizing multiplier size
        # alone tends to place the stationarity residual at the tolerance.
        bounds.append((0.0, tolerance))
        result = linprog(
            np.r_[np.zeros(inequalities.size), 1.0],
            A_ub=lhs if lhs.size else None,
            b_ub=rhs if rhs.size else None,
            bounds=bounds,
            method="highs",
            options={"primal_feasibility_tolerance": 1e-9, "dual_feasibility_tolerance": 1e-9},
        )
        lp_status = int(result.status)
        lp_message = str(result.message)
        if result.success:
            multipliers = result.x[:-1]
            source = "independent stationarity-margin LP"
        else:
            source = "zero diagnostic multipliers after LP did not return a witness"

    # Recompute all residuals using the actual design and returned multipliers.
    # The LP's feasibility status is insufficient to certify these tolerances.
    gradient_lagrangian = gradient + jacobian.T @ multipliers
    projected = design - np.clip(design - gradient_lagrangian, 0.0, 1.0)
    dual = float(np.max(-multipliers, initial=0.0))
    stationarity = float(np.max(np.abs(projected)))
    complementarity = float(np.max(np.abs(multipliers * inequalities), initial=0.0))
    verified = (
        primal <= primal_tolerance
        and dual <= dual_tolerance
        and stationarity <= tolerance
        and complementarity <= complementarity_tolerance
    )
    return KKTReport(
        verified=verified,
        multipliers=multipliers,
        multiplier_source=source,
        lp_status=lp_status,
        lp_message=lp_message,
        inequality_violation=inequality_violation,
        box_violation=box_violation,
        primal_violation=primal,
        dual_violation=dual,
        stationarity_inf=stationarity,
        complementarity_inf=complementarity,
        projected_residual=projected,
    )
