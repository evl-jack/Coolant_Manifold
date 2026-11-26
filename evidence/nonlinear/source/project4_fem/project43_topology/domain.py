"""Conforming planar FEM mesh for the sharp junction and its design region."""

from dataclasses import dataclass

import numpy as np

from ..fem import TaylorHood
from ..junction import boundary_data, strip_triangles
from .density import DesignGrid


REFERENCE_HEIGHT_MM = 4.0


@dataclass(frozen=True)
class TopologyDomain:
    fem: TaylorHood
    fixed: np.ndarray
    prescribed: np.ndarray
    ports: tuple[np.ndarray, np.ndarray, np.ndarray]
    flux: np.ndarray
    cell_index: np.ndarray
    max_spacing_mm: float


def _refined_axis(cuts, spacing_mm):
    cuts = np.unique(cuts)
    segments = []
    for lower, upper in zip(cuts[:-1], cuts[1:]):
        count = int(np.ceil((upper - lower) / spacing_mm))
        segments.append(np.linspace(lower, upper, count + 1)[:-1])
    return np.concatenate((*segments, cuts[-1:]))


def build_domain(grid: DesignGrid, max_spacing_mm: float) -> TopologyDomain:
    """Build a mesh of the original T union the design box.

    Every design-cell face follows FEM edges. The row-major cell index is
    constant over each element and is -1 in the fixed fluid outside the box.
    The spacing bound applies to the Cartesian axes before triangulation.
    """
    if not np.isfinite(max_spacing_mm) or max_spacing_mm <= 0.0:
        raise ValueError("The maximum FEM axis spacing must be finite and positive")

    x_axis = _refined_axis(np.r_[0.0, grid.x_edges_mm, 9.5, 10.5, 20.0], max_spacing_mm)
    y_axis = _refined_axis(np.r_[grid.y_edges_mm, 4.0, 7.0], max_spacing_mm)
    xx, yy = np.meshgrid(x_axis, y_axis)
    points_mm = np.column_stack((xx.ravel(), yy.ravel()))
    triangles = strip_triangles(np.arange(xx.size).reshape(xx.shape))

    x, y = points_mm[triangles].mean(axis=1).T
    original_fluid = (y < 4.0) | ((x > 9.5) & (x < 10.5))
    designable = (x > 6.0) & (x < 16.0) & (y < 6.0)
    keep = original_fluid | designable
    used_vertices, local_vertices = np.unique(triangles[keep], return_inverse=True)
    fem = TaylorHood(points_mm[used_vertices] / REFERENCE_HEIGHT_MM, local_vertices.reshape(-1, 3))

    cell_index = np.full(len(fem.triangles), -1, dtype=np.intp)
    inside = designable[keep]
    ix = np.searchsorted(grid.x_edges_mm, x[designable], side="right") - 1
    iy = np.searchsorted(grid.y_edges_mm, y[designable], side="right") - 1
    cell_index[inside] = np.ravel_multi_index((iy, ix), grid.shape)

    fixed, prescribed, ports, flux = boundary_data(fem)
    _validate_mesh(grid, fem, cell_index, original_fluid[keep], ports)
    return TopologyDomain(fem, fixed, prescribed, ports, flux, cell_index, max_spacing_mm)


def _validate_mesh(grid, fem, cell_index, original_fluid, ports):
    area = fem.weights.sum(axis=1) * REFERENCE_HEIGHT_MM**2
    designable = cell_index >= 0
    fixed_area = area[~designable].sum()
    cell_areas = np.bincount(cell_index[designable], weights=area[designable], minlength=grid.size)
    np.testing.assert_allclose(area.sum(), 101.0, rtol=1e-12)
    np.testing.assert_allclose(fixed_area, 41.0, rtol=1e-12)
    np.testing.assert_allclose(cell_areas, grid.cell_area_mm2, rtol=1e-12)
    np.testing.assert_allclose(area[original_fluid].sum(), 83.0, rtol=1e-12)
    np.testing.assert_array_equal(grid.base_raw[cell_index[designable]], original_fluid[designable])

    # No triangle may straddle the boundary of its assigned design cell.
    iy, ix = np.unravel_index(cell_index[designable], grid.shape)
    corners = fem.points[fem.triangles[designable]] * REFERENCE_HEIGHT_MM
    lower = np.column_stack((grid.x_edges_mm[ix], grid.y_edges_mm[iy]))
    upper = np.column_stack((grid.x_edges_mm[ix + 1], grid.y_edges_mm[iy + 1]))
    tolerance_mm = 1e-12
    if np.any(corners.min(axis=1) < lower - tolerance_mm) or np.any(corners.max(axis=1) > upper + tolerance_mm):
        raise ValueError("A FEM triangle crosses an assigned design-cell face")

    for edges, normal_axis, location, lower, upper in zip(
        ports, (0, 0, 1), (0.0, 20.0, 7.0), (0.0, 0.0, 9.5), (4.0, 4.0, 10.5)
    ):
        endpoints = fem.points[edges[:, :2]] * REFERENCE_HEIGHT_MM
        np.testing.assert_allclose(endpoints[..., normal_axis], location, rtol=0.0, atol=tolerance_mm)
        tangents = endpoints[..., 1 - normal_axis]
        np.testing.assert_allclose((tangents.min(), tangents.max()), (lower, upper), rtol=0.0, atol=tolerance_mm)
        widths = np.abs(tangents[:, 1] - tangents[:, 0])
        np.testing.assert_allclose(widths.sum(), upper - lower, rtol=1e-12)
