"""Tangent and adjoint algebra for a converged discrete residual R(z,m)=0.

The caller supplies the exact state Jacobian and partial derivatives after
boundary enforcement. Geometry and boundary lifting derivatives belong in
those partials. This module cannot manufacture them from a forward solve.

Sparse direct solves are an initial backend for small systems. Keep this
interface when replacing the backend with scalable preconditioned solvers.
"""

import warnings
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import MatrixRankWarning, spsolve


def _linear_solve(matrix, rhs):
    matrix = csc_matrix(matrix, dtype=float)
    rhs = np.asarray(rhs, dtype=float)
    if (matrix.shape[0] == 0 or matrix.shape[0] != matrix.shape[1]
            or rhs.shape != (matrix.shape[0],)):
        raise ValueError("Supply a nonempty square Jacobian and matching vector.")
    if not np.isfinite(matrix.data).all() or not np.isfinite(rhs).all():
        raise ValueError("Sensitivity equations must contain finite values.")
    with warnings.catch_warnings():
        warnings.simplefilter("error", MatrixRankWarning)
        solution = np.asarray(spsolve(matrix, rhs), dtype=float)
    if not np.isfinite(solution).all():
        raise ArithmeticError("The sensitivity solve returned nonfinite values.")
    return solution


def solve_tangent(state_jacobian, residual_parameter_direction):
    """Solve R_z * dz = -R_m*d for one supplied parameter direction d."""
    return _linear_solve(
        state_jacobian, -np.asarray(residual_parameter_direction, dtype=float)
    )


def solve_adjoint(state_jacobian, response_state_partial):
    """Solve R_z.T * lambda = J_z for one scalar response."""
    matrix = csc_matrix(state_jacobian, dtype=float)
    return _linear_solve(matrix.T, response_state_partial)


def adjoint_directional_derivative(
    response_parameter_direction, adjoint, residual_parameter_direction
):
    """Return J_m*d - lambda.T*(R_m*d), with free state held in the partials.

    Call this for each parameter direction or replace with a local assembled
    vector-Jacobian product for a large field. Do not form dense global
    state-by-design matrices solely to use this helper.
    """
    direct = float(response_parameter_direction)
    adjoint = np.asarray(adjoint, dtype=float)
    partial = np.asarray(residual_parameter_direction, dtype=float)
    if adjoint.ndim != 1 or partial.shape != adjoint.shape:
        raise ValueError("Adjoint and residual partial must be matching vectors.")
    if (not np.isfinite(direct) or not np.isfinite(adjoint).all()
            or not np.isfinite(partial).all()):
        raise ValueError("Derivative inputs must be finite.")
    result = float(direct - adjoint @ partial)
    if not np.isfinite(result):
        raise ArithmeticError("Nonfinite total directional derivative.")
    return result
