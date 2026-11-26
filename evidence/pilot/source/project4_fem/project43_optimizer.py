"""Scaled two-control SLSQP interface for an externally supplied flow response."""

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize


BOUNDS = ((0.0375, 1.0), (0.0, 1.0))


def to_physical(x):
    x = np.asarray(x, dtype=float)
    if x.shape != (2,) or not np.all(np.isfinite(x)):
        raise ValueError("Expected two finite normalized controls")
    return np.array([4.0 * x[0], x[1]])


def to_normalized(design_mm):
    design_mm = np.asarray(design_mm, dtype=float)
    if design_mm.shape != (2,) or not np.all(np.isfinite(design_mm)):
        raise ValueError("Expected two finite physical controls")
    return np.array([design_mm[0] / 4.0, design_mm[1]])


@dataclass(frozen=True)
class Evaluation:
    objective: float
    gradient: np.ndarray
    constraints: np.ndarray
    jacobian: np.ndarray
    power: float
    side_share: float


class DesignEvaluator:
    """response(h_mm, s_mm) returns power, share and physical derivatives.

    The response mapping requires power, side_share, power_gradient_mm,
    side_share_gradient_mm and converged. A failed response raises ValueError.
    """

    def __init__(self, response, reference_power):
        if not np.isfinite(reference_power) or reference_power <= 0:
            raise ValueError("Reference power must be positive and finite")
        self.response = response
        self.reference_power = float(reference_power)
        self.cache = {}

    @property
    def unique_evaluations(self):
        return len(self.cache)

    def evaluate(self, x):
        physical = to_physical(x)
        key = tuple(np.asarray(x, dtype=float))
        if key not in self.cache:
            data = self.response(*physical)
            if not data or data.get("converged") is not True:
                raise ValueError("Flow response did not converge")
            power = float(data["power"])
            share = float(data["side_share"])
            power_grad = np.asarray(data["power_gradient_mm"], dtype=float)
            share_grad = np.asarray(data["side_share_gradient_mm"], dtype=float)
            if (power_grad.shape != (2,) or share_grad.shape != (2,) or
                    not np.all(np.isfinite([power, share])) or
                    not np.all(np.isfinite(power_grad)) or
                    not np.all(np.isfinite(share_grad)) or power <= 0):
                raise ValueError("Invalid flow response or derivatives")
            chain = np.array([4.0, 1.0])
            ds = share_grad * chain / 0.001
            self.cache[key] = Evaluation(
                power / self.reference_power, power_grad * chain / self.reference_power,
                np.array([(share - .499) / .001, (.501 - share) / .001,
                          1 - key[0] - key[1]]),
                np.array([ds, -ds, [-1., -1.]]), power, share)
        return self.cache[key]


def optimize(evaluator, initial=(.3, .3)):
    """Return the SciPy result; the caller must inspect independent KKT data."""
    def objective(x):
        value = evaluator.evaluate(x)
        return value.objective, value.gradient

    return minimize(objective, np.asarray(initial, dtype=float), jac=True,
                    method="SLSQP", bounds=BOUNDS,
                    constraints=[{"type": "ineq",
                                  "fun": lambda x: evaluator.evaluate(x).constraints,
                                  "jac": lambda x: evaluator.evaluate(x).jacobian}],
                    options={"ftol": 1e-10, "maxiter": 35})
