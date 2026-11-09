"""Analytic partial derivatives of terms already specified in Analysis.md."""

import math


def brinkman_coefficient_and_partials(gamma, alpha_max, q_alpha):
    """Return alpha and partial derivatives at one supplied field value.

    Treat gamma, alpha_max, and q_alpha as independent for these partials.
    The caller must add chain-rule terms if these inputs are dependent.
    """
    if not math.isfinite(gamma) or not 0 <= gamma <= 1:
        raise ValueError("The design field must be finite and lie in [0, 1].")
    if not math.isfinite(alpha_max) or alpha_max < 0:
        raise ValueError("Maximum drag must be finite and nonnegative.")
    if not math.isfinite(q_alpha) or q_alpha <= 0:
        raise ValueError("The interpolation parameter must be finite and positive.")
    denominator = 1.0 + q_alpha * gamma
    result = {
        "alpha": alpha_max * (1.0 - gamma) / denominator,
        "d_alpha_d_gamma": -alpha_max * (1.0 + q_alpha) / denominator**2,
        "d_alpha_d_alpha_max": (1.0 - gamma) / denominator,
        "d_alpha_d_q_alpha": -alpha_max * gamma * (1.0 - gamma) / denominator**2,
    }
    if not all(math.isfinite(value) for value in result.values()):
        raise ArithmeticError("Nonfinite Brinkman coefficient or derivative.")
    return result


def plate_drag_and_partials(viscosity_pa_s, depth_m):
    """Depth-averaged laminar plate drag only; do not add this in resolved 3D."""
    if not math.isfinite(viscosity_pa_s) or viscosity_pa_s <= 0:
        raise ValueError("Dynamic viscosity must be finite and positive.")
    if not math.isfinite(depth_m) or depth_m <= 0:
        raise ValueError("Depth must be finite and positive.")
    result = {
        "coefficient": 12.0 * viscosity_pa_s / depth_m**2,
        "d_coefficient_d_viscosity": 12.0 / depth_m**2,
        "d_coefficient_d_depth": -24.0 * viscosity_pa_s / depth_m**3,
    }
    if not all(math.isfinite(value) for value in result.values()):
        raise ArithmeticError("Nonfinite plate-drag coefficient or derivative.")
    return result
