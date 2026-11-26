"""Independent Project 4.3 verification cases. Run as a module on Sol."""
import json
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from .project43_baseline import Geometry, TaylorHood, baseline, boundary_data, fitted_mesh
from .verify import channel


def geometry():
    g = Geometry()
    fem = TaylorHood(*fitted_mesh(g, 0.25))
    fixed, values, ports, flux = boundary_data(fem)
    points = fem.points * g.height_mm
    edges = np.concatenate([fem.triangles[:, [0, 1]], fem.triangles[:, [1, 2]], fem.triangles[:, [2, 0]]])
    graph = coo_matrix((np.ones(len(edges)), edges.T), shape=(len(points), len(points)))
    components = connected_components(graph, directed=False, return_labels=False)
    assert components == 1
    widths = [np.linalg.norm(points[p[:, 0]] - points[p[:, 1]], axis=1).sum() for p in ports]
    area = fem.weights.sum() * g.height_mm**2
    np.testing.assert_allclose([area, *widths], [83, 4, 4, 1], rtol=1e-12)
    midpoints = fem.velocity_points[fem.boundary_edges[:, 2]] * g.height_mm
    x, y = midpoints.T
    outlet = np.isclose(x, 20) | np.isclose(y, 7)
    constrained = np.unique(fem.boundary_edges[~outlet])
    np.testing.assert_array_equal(np.sort(fixed), np.sort(np.r_[constrained, fem.n_u + constrained]))
    state = np.zeros(fem.size)
    state[fixed] = values
    x, y = fem.velocity_points[constrained].T
    np.testing.assert_allclose(state[constrained], np.where(np.isclose(x, 0), 6*y*(1-y), 0))
    np.testing.assert_allclose(state[constrained + fem.n_u], 0)
    np.testing.assert_allclose(flux[0] @ state, 1., rtol=1e-12)
    return dict(area_mm2=area, port_widths_mm=widths, connected_components=int(components))


def laminar():
    result, *_ = baseline(Geometry(), 0.25, 50.)
    reference = channel(16, 16, inertia=50.)[-1]
    assert reference["velocity_l2"] < 1e-8 and reference["relative_flux_error"] < 1e-8
    return dict(junction=result, clear_channel=reference)


def brinkman():
    result = channel(16, 16, alpha=25.)[-1]
    assert result["velocity_l2"] < 1e-3 and result["relative_flux_error"] < 1e-4
    assert result["relative_mass_imbalance"] < 1e-10
    return result


if __name__ == "__main__":
    case = sys.argv[1]
    result = {"geometry": geometry, "laminar": laminar, "brinkman": brinkman}[case]()
    print(json.dumps(dict(case=case, result=result), indent=2, allow_nan=False))
