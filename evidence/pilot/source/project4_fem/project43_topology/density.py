"""Common cell-centered design density and reverse sensitivities."""

from dataclasses import dataclass

import numpy as np
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import splu


DESIGN_X_MM = (6., 16.)
DESIGN_Y_MM = (0., 6.)
FIXED_FLUID_AREA_MM2 = 41.
REDUCED_SIZE = 20
_BOX_ROUNDOFF = 32 * np.finfo(float).eps


def _unit_vector(values, size, name):
    values = np.asarray(values, dtype=float)
    if (values.shape != (size,) or not np.all(np.isfinite(values))
            or np.any(values < -_BOX_ROUNDOFF) or np.any(values > 1 + _BOX_ROUNDOFF)):
        raise ValueError(f"{name} must be a finite length-{size} vector in [0, 1]")
    return values


class DesignGrid:
    """Uniform design cells. Vectors use y-outer, x-inner ordering."""

    def __init__(self, spacing_mm):
        if not np.isfinite(spacing_mm) or spacing_mm <= 0:
            raise ValueError("Design spacing must be finite and positive")
        # The baseline interfaces and exterior of the design region must be,
        # by construction, cell faces at every supported design spacing.
        offsets = np.array([10., 6., 3.5, 4.5, 4.])
        counts = offsets / spacing_mm
        if not np.allclose(counts, np.rint(counts), rtol=0, atol=1e-10):
            raise ValueError("Design spacing must align with the sharp T-junction and design bounds")
        nx, ny, left, right, roof = (int(count) for count in np.rint(counts))
        self.spacing_mm = float(spacing_mm)
        self.shape = (ny, nx)
        self.size = ny * nx
        self.cell_area_mm2 = self.spacing_mm ** 2
        self.x_edges_mm = np.linspace(*DESIGN_X_MM, nx + 1)
        self.y_edges_mm = np.linspace(*DESIGN_Y_MM, ny + 1)
        x_centers = (self.x_edges_mm[1:] + self.x_edges_mm[:-1]) / 2
        y_centers = (self.y_edges_mm[1:] + self.y_edges_mm[:-1]) / 2
        x, y = np.meshgrid(x_centers, y_centers, indexing="xy")
        self.centers_mm = np.column_stack((x.ravel(), y.ravel()))
        base = np.zeros(self.shape)
        base[:roof, :] = 1.
        base[roof:, left:right] = 1.
        self.base_raw = base.ravel()


def _wendland_basis(points, y_centers, radius_mm):
    x, y = np.meshgrid([7., 9., 11., 13., 15.], y_centers, indexing="xy")
    centers = np.column_stack((x.ravel(), y.ravel()))
    r = np.linalg.norm(points[:, None, :] - centers, axis=-1) / radius_mm
    values = np.maximum(0., 1. - r) ** 4 * (4. * r + 1.)
    return values / np.maximum(1., values.sum(axis=1, keepdims=True))


class ReducedDensity:
    """Ten fluid expansions followed by ten fluid contractions."""

    def __init__(self, grid):
        self.grid = grid
        expansion = _wendland_basis(grid.centers_mm, [4.5, 5.5], 2.)
        contraction = _wendland_basis(grid.centers_mm, [1., 3.], 2.5)
        base = grid.base_raw[:, None]
        self.jacobian = np.hstack(((1. - base) * expansion, -base * contraction))

    def raw(self, controls):
        controls = _unit_vector(controls, REDUCED_SIZE, "Reduced controls")
        return self.grid.base_raw + self.jacobian @ controls

    def pullback(self, raw_gradient):
        return self.jacobian.T @ raw_gradient


def _neumann_stiffness(n):
    # On square cells each internal face has length / center distance = 1.
    degree = np.full(n, 2.)
    degree[[0, -1]] = 1.
    return diags((-np.ones(n - 1), degree, -np.ones(n - 1)), [-1, 0, 1], format="csc")


def project(filtered, beta, threshold):
    if beta == 0:
        return filtered.copy(), np.ones_like(filtered)
    left = np.tanh(beta * threshold)
    denominator = left + np.tanh(beta * (1. - threshold))
    transition = beta * (filtered - threshold)
    projected = (left + np.tanh(transition)) / denominator
    tail = np.exp(-2. * np.abs(transition))
    derivative = beta * (4. * tail / (1. + tail) ** 2) / denominator
    return projected, derivative


def brinkman(projected, q, alpha_max):
    """Dimensional resistance and derivative with respect to projected density."""
    denominator = q + projected
    resistance = alpha_max * q * (1. - projected) / denominator
    derivative = -alpha_max * q * (q + 1.) / denominator ** 2
    return resistance, derivative


@dataclass(frozen=True)
class DensityState:
    raw: np.ndarray
    filtered: np.ndarray
    projected: np.ndarray
    resistance: np.ndarray
    projection_derivative: np.ndarray
    resistance_derivative: np.ndarray
    area_mm2: float


class DensityPipeline:
    """One factorization shared by forward and reverse Helmholtz solves."""

    def __init__(self, grid, filter_radius_mm=0.5, threshold=0.5):
        if not np.isfinite(filter_radius_mm) or filter_radius_mm <= 0:
            raise ValueError("Filter radius must be finite and positive")
        if not 0 < threshold < 1:
            raise ValueError("Projection threshold must be in (0, 1)")
        self.grid = grid
        self.filter_radius_mm = float(filter_radius_mm)
        self.threshold = float(threshold)
        ny, nx = grid.shape
        stiffness = (kron(eye(ny), _neumann_stiffness(nx))
                     + kron(_neumann_stiffness(ny), eye(nx)))
        mass = grid.cell_area_mm2 * eye(grid.size, format="csc")
        helmholtz = mass + self.filter_radius_mm ** 2 * stiffness
        self._factor = splu(helmholtz.tocsc())

    def filter(self, raw):
        return self._factor.solve(self.grid.cell_area_mm2 * raw)

    def filter_pullback(self, filtered_gradient):
        return self.grid.cell_area_mm2 * self._factor.solve(filtered_gradient, trans="T")

    def evaluate(self, raw, q, beta, alpha_max):
        raw = _unit_vector(raw, self.grid.size, "Raw density").copy()
        if not np.isfinite([q, beta, alpha_max]).all() or q <= 0 or beta < 0 or alpha_max <= 0:
            raise ValueError("q and alpha_max must be positive and beta nonnegative, all finite")
        filtered = self.filter(raw)
        projected, projection_derivative = project(filtered, beta, self.threshold)
        resistance, resistance_derivative = brinkman(projected, q, alpha_max)
        area_mm2 = FIXED_FLUID_AREA_MM2 + self.grid.cell_area_mm2 * projected.sum()
        return DensityState(raw, filtered, projected, resistance, projection_derivative,
                            resistance_derivative, float(area_mm2))

    def pullback(self, state, resistance_gradient):
        filtered_gradient = state.projection_derivative * state.resistance_derivative * resistance_gradient
        return self.filter_pullback(filtered_gradient)

    def area_gradient(self, state):
        return self.filter_pullback(self.grid.cell_area_mm2 * state.projection_derivative)
