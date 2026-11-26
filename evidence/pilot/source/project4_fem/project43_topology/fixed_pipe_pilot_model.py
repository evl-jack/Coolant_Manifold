"""Fixed-pipe Brinkman pilot built on the existing topology FEM and adjoints."""

import numpy as np
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import splu
from shapely.geometry import box

from ..fem import TaylorHood
from ..junction import boundary_data, strip_triangles
from .density import DesignGrid, _neumann_stiffness, _unit_vector
from .domain import REFERENCE_HEIGHT_MM, TopologyDomain, _refined_axis
from .evaluator import HEIGHT_M, REYNOLDS, VELOCITY_SCALE, Evaluator
from .reconstruction import Reconstruction, reconstruct
from .wall_displacement_density_v1 import WallDisplacementDensityV1
from .wall_displacement_geometry_v1 import AREA_ROUNDOFF_MM2, _clipped, _turning_report, sharp_geometry


DESIGN_BOUNDS_MM = (8., 4., 12., 6.)
MAIN_BOUNDS_MM = (0., 0., 20., 4.)
COLLAR_BOUNDS_MM = (9.5, 6., 10.5, 7.)
SPACING_MM = .25
FILTER_RADIUS_MM = .5
SIDE_CONSTRAINT_SCALE = .001
EDITABLE_COUNT = 128


def _in_box(x, y, bounds):
    xmin, ymin, xmax, ymax = bounds
    return (xmin < x) & (x < xmax) & (ymin < y) & (y < ymax)


class FixedPipeDensity(WallDisplacementDensityV1):
    """Keep the established filter solve and transpose with fixed pipe anchors."""

    def __init__(self, grid):
        self.grid = grid
        self.filter_radius_mm = FILTER_RADIUS_MM
        self.threshold = .5
        self.free = _in_box(*grid.centers_mm.T, DESIGN_BOUNDS_MM)
        self.free.setflags(write=False)
        self.size = int(self.free.sum())
        if self.size != EDITABLE_COUNT:
            raise ValueError(f"Expected {EDITABLE_COUNT} editable cells, got {self.size}")
        self._base_raw = grid.base_raw.copy()

        ny, nx = grid.shape
        radius_squared = self.filter_radius_mm ** 2
        stiffness = kron(eye(ny), _neumann_stiffness(nx)) + kron(_neumann_stiffness(ny), eye(nx))
        mass = grid.cell_area_mm2 * eye(grid.size, format="csc")
        helmholtz = (mass + radius_squared * stiffness).tocsr()
        rows = helmholtz[self.free]
        self._fixed_rhs = -(rows[:, ~self.free] @ self._base_raw[~self.free])

        # The standard grid ends at y=6. Add the face against the fixed
        # collar/solid cells immediately above it, as an eliminated halo row.
        x, y = grid.centers_mm[self.free].T
        top = y == grid.y_edges_mm[-1] - grid.spacing_mm / 2.
        collar = _in_box(x, y + grid.spacing_mm, COLLAR_BOUNDS_MM)
        self._fixed_rhs += radius_squared * (top & collar)
        self._factor = splu((rows[:, self.free] + diags(radius_squared * top)).tocsc())

    def evaluate(self, raw, q, beta, alpha_max):
        density = super().evaluate(raw, q, beta, alpha_max)
        fixed = ~self.free
        for name in ("raw", "filtered", "projected"):
            values = _unit_vector(getattr(density, name), self.grid.size, name)
            if not np.array_equal(values[fixed], self._base_raw[fixed]):
                raise RuntimeError(f"{name} density changed a fixed cell")
        if not np.isfinite(density.resistance).all() or np.any(density.resistance[fixed & (self._base_raw == 1.)]):
            raise RuntimeError("A fixed fluid cell acquired Brinkman resistance or the field is nonfinite")
        return density


def build_fixed_pipe():
    """Use the common structured mesh machinery on the exact pilot envelope."""
    grid = DesignGrid(SPACING_MM)
    pipeline = FixedPipeDensity(grid)
    x_axis = _refined_axis(np.r_[0., grid.x_edges_mm, 20.], SPACING_MM)
    y_axis = _refined_axis(np.r_[grid.y_edges_mm, 7.], SPACING_MM)
    xx, yy = np.meshgrid(x_axis, y_axis)
    points = np.column_stack((xx.ravel(), yy.ravel()))
    triangles = strip_triangles(np.arange(xx.size).reshape(xx.shape))
    x, y = points[triangles].mean(axis=1).T
    editable = _in_box(x, y, DESIGN_BOUNDS_MM)
    keep = editable | _in_box(x, y, MAIN_BOUNDS_MM) | _in_box(x, y, COLLAR_BOUNDS_MM)
    used, compact = np.unique(triangles[keep], return_inverse=True)
    fem = TaylorHood(points[used] / REFERENCE_HEIGHT_MM, compact.reshape(-1, 3))
    fem.inertia = REYNOLDS
    fixed, prescribed, ports, flux = boundary_data(fem)

    cell_index = np.full(len(fem.triangles), -1, dtype=np.intp)
    ix = np.searchsorted(grid.x_edges_mm, x[editable], side="right") - 1
    iy = np.searchsorted(grid.y_edges_mm, y[editable], side="right") - 1
    cell_index[editable[keep]] = np.ravel_multi_index((iy, ix), grid.shape)

    area = fem.weights.sum(axis=1) * REFERENCE_HEIGHT_MM ** 2
    in_design = cell_index >= 0
    cell_areas = np.bincount(cell_index[in_design], weights=area[in_design], minlength=grid.size)
    np.testing.assert_allclose(cell_areas, grid.cell_area_mm2 * pipeline.free, rtol=1e-12)
    np.testing.assert_allclose(area[~in_design].sum(), 81., rtol=1e-12)
    lengths = [np.linalg.norm(fem.points[edges[:, 1]] - fem.points[edges[:, 0]], axis=1).sum()
               * REFERENCE_HEIGHT_MM for edges in ports]
    np.testing.assert_allclose(lengths, [4., 4., 1.], rtol=1e-12)
    return TopologyDomain(fem, fixed, prescribed, ports, flux, cell_index, SPACING_MM), pipeline


class FixedPipeEvaluator(Evaluator):
    """The original nonlinear state and adjoints with one reference-flow inequality."""

    def __init__(self, domain, pipeline, stage, power_reference, side_reference, *, initial_state=None):
        super().__init__(domain, pipeline, stage, power_reference, initial_state=initial_state)
        if not np.isfinite(side_reference):
            raise ValueError("Side-flow reference must be finite")
        self.side_reference = float(side_reference)

    def design_constraints(self, density, fraction, fraction_gradient):
        return (np.array([(self.side_reference - fraction) / SIDE_CONSTRAINT_SCALE]),
                -fraction_gradient[None, :] / SIDE_CONSTRAINT_SCALE)

    def solve(self, resistance, initial):
        state = super().solve(resistance, initial)
        residual = state.residual_history[-1]
        if not np.isfinite(residual) or residual > 1e-10 or not np.isfinite(state.values).all():
            raise RuntimeError("The nonlinear state failed basic finiteness or convergence checks")
        return state

    def evaluate(self, design):
        result = super().evaluate(design)
        if (not np.isfinite(result.objective) or result.objective <= 0.
                or not np.isfinite(result.gradient).all() or not np.isfinite(result.jacobian).all()):
            self.last_design = None
            raise RuntimeError("The pumping power or adjoint gradients failed basic validity checks")
        return result


def original_geometry():
    geometry = sharp_geometry()
    return Reconstruction(geometry, (geometry,), (0, 0, 0), ())


def reconstruct_geometry(pipeline, density):
    """Clip the threshold contour to the transition and retain every fixed wall."""
    original = sharp_geometry()
    rectangle = box(*DESIGN_BOUNDS_MM)
    result = _clipped(original, rectangle, reconstruct(pipeline.grid, density.projected))
    outside_change = result.geometry.symmetric_difference(original).difference(rectangle).area
    if (outside_change > AREA_ROUNDOFF_MM2 or not result.geometry.covers(box(*MAIN_BOUNDS_MM))
            or not result.geometry.covers(box(*COLLAR_BOUNDS_MM))):
        raise RuntimeError("Reconstruction changed the fixed main pipe or outlet collar")
    return result


def metrics_for(evaluation, geometry, pipeline):
    """Combine immersed flow scalars with geometric contour diagnostics."""
    turning = _turning_report(sharp_geometry(), geometry, box(*DESIGN_BOUNDS_MM))
    inlet, main, side = evaluation.flux
    inlet_flow = inlet * HEIGHT_M * VELOCITY_SCALE
    return dict(power_w_per_m=evaluation.power_w_per_m, normalized_power=evaluation.objective,
                hydraulic_resistance_pa_s_per_m2=float(evaluation.power_w_per_m / inlet_flow ** 2),
                inlet_flow_m2_per_s=float(inlet_flow), main_fraction=float(main / inlet),
                side_fraction=float(side / inlet), physical_area_mm2=geometry.area_mm2,
                max_turn_degrees=turning["max_turn_degrees"], turn_secant_mm=turning["secant_length_mm"],
                turn_sample_spacing_mm=turning["arc_step_mm"], ports_connected=geometry.port_component is not None,
                editable_cells=pipeline.size, mass_imbalance=evaluation.mass_imbalance,
                nonlinear_residual=evaluation.nonlinear_residual)
