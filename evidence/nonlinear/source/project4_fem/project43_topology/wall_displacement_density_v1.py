"""Compact wall-displacement controls for the original sharp T-junction."""

import numpy as np
from scipy.sparse import eye, kron
from scipy.sparse.linalg import splu

from .density import DESIGN_X_MM, DESIGN_Y_MM, DensityPipeline, _neumann_stiffness, _unit_vector


WALL_SEGMENTS_MM = (
    ((0., 0.), (20., 0.)),
    ((0., 4.), (9.5, 4.)),
    ((10.5, 4.), (20., 4.)),
    ((9.5, 4.), (9.5, 7.)),
    ((10.5, 4.), (10.5, 7.)),
)
EDITABLE_BOX_MM = (6.5, 0., 15.5, 5.5)
FLUID_CORE_BOXES_MM = ((0., 1., 20., 3.), (9.5, 2., 10.5, 7.))
WALL_DISTANCE_MM = 1.


def wall_displacement_mask(grid, core_boxes=FLUID_CORE_BOXES_MM):
    """Identify cells wholly eligible under the v1 wall and core rules."""
    x, y = grid.x_edges_mm, grid.y_edges_mm
    if (x[0], x[-1]) != DESIGN_X_MM or (y[0], y[-1]) != DESIGN_Y_MM:
        raise ValueError("Wall displacement v1 requires the original T-junction design envelope")
    x0, y0 = np.meshgrid(x[:-1], y[:-1], indexing="xy")
    x1, y1 = np.meshgrid(x[1:], y[1:], indexing="xy")
    xmin, ymin, xmax, ymax = EDITABLE_BOX_MM
    in_box = (x0 >= xmin) & (x1 <= xmax) & (y0 >= ymin) & (y1 <= ymax)
    outside_cores = np.ones(grid.shape, dtype=bool)
    for xmin, ymin, xmax, ymax in core_boxes:
        # Boundary contact has zero area and does not intersect the cell interior.
        overlap = (x0 < xmax) & (x1 > xmin) & (y0 < ymax) & (y1 > ymin)
        outside_cores &= ~overlap

    corners = np.stack((np.stack((x0, y0), axis=-1), np.stack((x1, y0), axis=-1),
                        np.stack((x1, y1), axis=-1), np.stack((x0, y1), axis=-1)), axis=-2)
    near_wall = np.zeros(grid.shape, dtype=bool)
    for start, end in np.asarray(WALL_SEGMENTS_MM):
        segment = end - start
        offset = corners - start
        parameter = np.clip((offset @ segment) / (segment @ segment), 0., 1.)
        separation = offset - parameter[..., None] * segment
        corner_distance_squared = np.sum(separation**2, axis=-1)
        near_wall |= np.all(corner_distance_squared <= WALL_DISTANCE_MM**2, axis=-1)
    return (in_box & outside_cores & near_wall).ravel()


class WallDisplacementDensityV1(DensityPipeline):
    """Solve the original Helmholtz filter with fixed binary cell values.

    The editable and fixed index sets are E and F. Partition the original
    H = M + r_f^2 K to solve H_EE f_E = M_E m_E - H_EF m_F, with f_F = m_F
    fixed to the binary sharp-junction values. The transpose sensitivity for
    a full filtered-field cotangent g is M_E H_EE^{-T} g_E in compact order.
    """

    core_boxes = FLUID_CORE_BOXES_MM

    def __init__(self, grid, filter_radius_mm=.5, threshold=.5):
        if not np.isfinite(filter_radius_mm) or filter_radius_mm <= 0.:
            raise ValueError("Filter radius must be finite and positive")
        if threshold != .5:
            raise ValueError("Wall displacement v1 uses the fixed 0.5 threshold")
        self.grid = grid
        self.filter_radius_mm = float(filter_radius_mm)
        self.threshold = float(threshold)
        self.free = wall_displacement_mask(grid, self.core_boxes)
        self.free.setflags(write=False)
        self.size = int(np.count_nonzero(self.free))
        if self.size == 0:
            raise ValueError("The design grid contains no editable wall cells")
        self._base_raw = _unit_vector(grid.base_raw, grid.size, "Sharp T raw density").copy()
        if not np.all((self._base_raw == 0.) | (self._base_raw == 1.)):
            raise ValueError("The fixed anchors require an original binary T-junction design")

        ny, nx = grid.shape
        stiffness = (kron(eye(ny), _neumann_stiffness(nx))
                     + kron(_neumann_stiffness(ny), eye(nx)))
        mass = grid.cell_area_mm2 * eye(grid.size, format="csc")
        helmholtz = (mass + self.filter_radius_mm**2 * stiffness).tocsr()
        editable_rows = helmholtz[self.free]
        self._fixed_rhs = -(editable_rows[:, ~self.free] @ self._base_raw[~self.free])
        self._factor = splu(editable_rows[:, self.free].tocsc())

    @property
    def initial(self):
        """Return the compact binary sharp T initialization in grid order."""
        return self._base_raw[self.free]

    def raw(self, design):
        """Expand controls, or validate a full vector supplied for replay."""
        design = np.asarray(design, dtype=float)
        if design.shape == (self.grid.size,):
            full = _unit_vector(design, self.grid.size, "Raw density")
            if not np.array_equal(full[~self.free], self._base_raw[~self.free]):
                raise ValueError("An inactive wall-design cell differs from its binary anchor")
            return full.copy()
        editable = _unit_vector(design, self.size, "Editable wall density")
        full = self._base_raw.copy()
        full[self.free] = editable
        return full

    def evaluate(self, raw, q, beta, alpha_max):
        """Apply the common projection and Brinkman map to the anchored field."""
        return super().evaluate(self.raw(raw), q, beta, alpha_max)

    def filter(self, raw):
        full = self.raw(raw)
        rhs = self.grid.cell_area_mm2 * full[self.free] + self._fixed_rhs
        full[self.free] = self._factor.solve(rhs)
        return full

    def filter_pullback(self, filtered_gradient):
        filtered_gradient = np.asarray(filtered_gradient, dtype=float)
        if filtered_gradient.shape != (self.grid.size,) or not np.isfinite(filtered_gradient).all():
            raise ValueError("Expected one finite filtered-field cotangent per design cell")
        return self.grid.cell_area_mm2 * self._factor.solve(filtered_gradient[self.free], trans="T")
