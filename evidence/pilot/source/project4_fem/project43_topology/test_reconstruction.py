"""Sol fixtures for topology reconstruction and the fitted-wall mesh."""

import numpy as np
import pytest

from project4_fem.project43_topology.density import DesignGrid
from project4_fem.project43_topology.reconstruction import (
    _cell_regions, fitted_mesh, reconstruct, write_reconstruction,
)


_UNIT_SQUARE = [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]


@pytest.mark.parametrize("values, rings, segments", [
    ([.49, -.10, .03, -.10], 1, 2),
    ([.49, -.10, .015, -.10], 2, 2),
    ([-.10, .49, -.10, .03], 1, 2),
    ([-.10, .49, -.10, .015], 2, 2),
    ([.25, -.25, .25, -.25], 2, 4),
    ([-.25, .25, -.25, .25], 2, 4),
])
def test_bilinear_saddle(values, rings, segments):
    actual_rings, actual_segments = _cell_regions(_UNIT_SQUARE, np.array(values))
    assert len(actual_rings) == rings
    assert len(actual_segments) == segments
    if segments == 4:
        assert all(segment[1] == (.5, .5) for segment in actual_segments)
        assert set(actual_rings[0]) & set(actual_rings[1]) == {(.5, .5)}


def test_reconstruction_extrema():
    grid = DesignGrid(.5)
    solid = reconstruct(grid, np.zeros(grid.size))
    assert solid.area_mm2 == 41.
    assert not solid.contour_segments_mm
    assert solid.port_component is None
    assert len(set(solid.port_components)) == 3

    fluid = reconstruct(grid, np.ones(grid.size))
    assert fluid.area_mm2 == 101.
    assert not fluid.contour_segments_mm
    assert fluid.port_component is not None
    assert len(fluid.polygons) == 1


def test_exact_saddle_does_not_connect_fluid_components():
    grid = DesignGrid(.5)
    projected = np.zeros(grid.size)
    x, y = grid.centers_mm.T
    projected[((x == 12.25) & (y == 5.25)) | ((x == 12.75) & (y == 5.75))] = 1.
    result = reconstruct(grid, projected)
    pockets = [polygon for i, polygon in enumerate(result.polygons) if i not in result.port_components]
    assert len(pockets) == 2
    contact = pockets[0].intersection(pockets[1])
    assert contact.geom_type == "Point"
    assert tuple(contact.coords[0]) == (12.5, 5.5)


def test_reconstruction_preserves_island_and_pocket(tmp_path):
    grid = DesignGrid(.5)
    raw = grid.base_raw.copy()
    x, y = grid.centers_mm.T
    raw[(x == 12.25) & (y == 2.25)] = 0.
    raw[(x == 13.25) & (y == 5.25)] = 1.
    result = reconstruct(grid, raw)
    assert len(result.polygons) == 2
    assert result.port_component is not None
    assert len(result.port_component.interiors) == 1
    path = tmp_path / "geometry.json"
    report = write_reconstruction(result, path)
    assert path.exists()
    assert report["connected_ports"]
    assert len(report["geometry_sha256"]) == 64
    assert sum(component["solid_islands"] for component in report["components"]) == 1


def test_fitted_mesh_reuse_fixed_contour():
    grid = DesignGrid(.5)
    result = reconstruct(grid, grid.base_raw)
    assert result.port_component is not None
    counts = []
    for spacing in (.5, .25):
        points, triangles = fitted_mesh(result.port_component, spacing)
        counts.append(len(triangles))
        np.testing.assert_allclose(points.min(axis=0), (0., 0.))
        np.testing.assert_allclose(points.max(axis=0), (5., 1.75))
    assert counts[1] > counts[0] > 0
