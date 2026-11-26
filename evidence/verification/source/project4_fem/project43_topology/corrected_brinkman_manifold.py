"""Paired fixed-mask manifold check of the planar displacement correction.

Let S be the nominal solid and ell = sqrt(mu / alpha). The corrected binary
mask adds (offset(S, +ell) intersect Omega) minus S. Straight segments are
translated into nominal fluid by ell. Mitered corners are a specified geometric
extension of the planar correction and have no independent physical claim.

The nominal solid, fitted reference, ports and existing mesh are held fixed.
In particular, leakage is always evaluated across the same nominal interface.
The correction enters only the drag term, where clipped subcells give exact
degree-four integration of the P2 velocity shape-function products. The
verified Taylor–Hood and Newton implementations are used unchanged. This
fixed binary verification does not provide a density derivative or adopt a
correction for production optimization.
"""

import argparse
from copy import copy
import json
import os
from pathlib import Path
from time import perf_counter
import traceback

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix

from ..fem import quadratic_basis, triangle_rule
from .density import DesignGrid
from .domain import REFERENCE_HEIGHT_MM, build_domain
from .evaluator import DesignTransfer, RESISTANCE_SCALE, VISCOSITY
from .gate2 import _baseline_record, _continuation
from .gate3_screen import PENALTIES, obstruction
from .leakage import binary_leakage
from .local_refinement import build_local_binary_domain


def _polygon_parts(geometry):
    if geometry.geom_type == "Polygon":
        if not geometry.is_empty:
            yield geometry
    elif geometry.geom_type in ("MultiPolygon", "GeometryCollection"):
        for item in geometry.geoms:
            yield from _polygon_parts(item)


def _integrate_clipped_mass(element, polygon, barycentric_nodes, quadrature_weights):
    """Integrate on signed fan triangles in original FE coordinates.

    Fan triangles can overlap or have negative orientation for concave rings.
    Signed Jacobians recover the region integral, with hole rings subtracted.
    The degree-four integrand is the product of two parent P2 basis functions.
    """
    from shapely.geometry.polygon import orient

    anchor = element[0]
    inverse = np.linalg.inv((element[1:] - anchor).T)
    mass = np.zeros((6, 6))
    for part in _polygon_parts(polygon):
        oriented = orient(part, sign=1.)
        for ring in (oriented.exterior, *oriented.interiors):
            points = np.asarray(ring.coords)
            first, second = points[1:-2], points[2:-1]
            ring_anchor = points[0]
            sides1, sides2 = first - ring_anchor, second - ring_anchor
            determinant = sides1[:, 0] * sides2[:, 1] - sides1[:, 1] * sides2[:, 0]
            nodes = (ring_anchor + barycentric_nodes[:, 1, None] * sides1[:, None, :]
                     + barycentric_nodes[:, 2, None] * sides2[:, None, :])
            xy = (nodes - anchor) @ inverse.T
            bary = np.concatenate((1 - xy.sum(axis=-1, keepdims=True), xy), axis=-1)
            basis, _ = quadratic_basis(bary)
            mass += np.einsum("e,q,eqa,eqb->ab", determinant, quadrature_weights, basis, basis)
    return mass


def _nominal_geometry(grid, raw):
    import shapely

    scale = REFERENCE_HEIGHT_MM
    x, y = grid.x_edges_mm / scale, grid.y_edges_mm / scale
    iy, ix = np.unravel_index(np.flatnonzero(raw == 0.), grid.shape)
    solid = shapely.union_all(shapely.box(x[ix], y[iy], x[ix + 1], y[iy + 1]))
    box = shapely.box(x[0], y[0], x[-1], y[-1])
    sharp = shapely.box(0., 0., 20. / scale, 4. / scale).union(
        shapely.box(9.5 / scale, 4. / scale, 10.5 / scale, 7. / scale))
    computational_domain = box.union(sharp)
    return solid, computational_domain


def assemble_added_drag(domain, grid, raw, alpha, shift_layers=1.):
    """Return the matrix for the added strip, plus assembly diagnostics."""
    import shapely
    from shapely.geometry import LineString, Polygon

    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError("Require a finite positive dimensional resistance")
    if not np.isfinite(shift_layers) or not (0 <= shift_layers <= 1):
        raise ValueError("Require a shift between zero and one Brinkman length")
    raw = np.asarray(raw)
    if raw.shape != (grid.size,) or not np.all((raw == 0.) | (raw == 1.)):
        raise ValueError("Require one binary raw value per design cell")
    fem = domain.fem
    area_scale = REFERENCE_HEIGHT_MM**2
    solid, computational_domain = _nominal_geometry(grid, raw)
    nominal_solid = np.zeros(len(fem.triangles), dtype=bool)
    designable = domain.cell_index >= 0
    nominal_solid[designable] = raw[domain.cell_index[designable]] == 0.
    np.testing.assert_allclose(fem.weights[nominal_solid].sum(), solid.area, rtol=1e-12, atol=1e-13)
    np.testing.assert_allclose(fem.weights.sum(), computational_domain.area, rtol=1e-12, atol=1e-13)
    nominal_fluid_area = (computational_domain.area - solid.area) * area_scale
    np.testing.assert_allclose(nominal_fluid_area, 83., rtol=1e-12, atol=1e-12)

    layer_mm = 1000 * np.sqrt(VISCOSITY / alpha)
    shift_mm = shift_layers * layer_mm
    diagnostics = dict(brinkman_length_mm=float(layer_mm), shift_mm=float(shift_mm),
                       shift_layers=float(shift_layers), nominal_fluid_area_mm2=float(nominal_fluid_area),
                       nominal_solid_area_mm2=float(solid.area * area_scale), corners="mitered")
    if not shift_mm:
        return csr_matrix((fem.size, fem.size)), {
            **diagnostics, "added_penalty_area_mm2": 0., "integrated_added_penalty_area_mm2": 0.,
            "elements_with_added_drag": 0,
        }

    shifted = solid.buffer(shift_mm / REFERENCE_HEIGHT_MM, join_style="mitre", mitre_limit=2.)
    extension = shifted.intersection(computational_domain).difference(solid)
    for edges in domain.ports:
        for edge in edges:
            if not extension.disjoint(LineString(fem.points[edge[:2]])):
                raise ValueError("The displacement-corrected penalty mask reaches a fixed aperture")

    # Only original fluid elements can contribute additional drag. The nominal
    # binary mask continues to be assembled by TaylorHood.drag_matrix.
    unique = np.flatnonzero(~nominal_solid)
    vertices = fem.points[fem.triangles[unique]]
    triangles = shapely.polygons(vertices)
    intersects = shapely.intersects(triangles, extension)
    clipped = shapely.intersection(triangles[intersects], extension)
    positive = shapely.area(clipped) > 0.
    element_ids = unique[intersects][positive]
    poly_nodes, poly_weights = triangle_rule()
    masses = np.array([_integrate_clipped_mass(fem.points[fem.triangles[index]], region,
                                               poly_nodes, poly_weights)
                       for index, region in zip(element_ids, clipped[positive])])
    if not len(element_ids):
        raise RuntimeError("The positive displacement adds no integrated drag")
    integrated_area = float(masses.sum())
    np.testing.assert_allclose(integrated_area, extension.area, rtol=1e-10, atol=1e-13)
    np.testing.assert_allclose(masses.sum(axis=(1, 2)), shapely.area(clipped[positive]),
                               rtol=1e-9, atol=1e-13)

    # A whole-triangle fit must reproduce the stock P2 drag integral. Use the
    # same clipped-element geometry to exercise exactly the production basis.
    i = element_ids[0]
    element = fem.points[fem.triangles[i]]
    full = _integrate_clipped_mass(element, Polygon(element), poly_nodes, poly_weights)
    stock = np.einsum("q,qa,qb->ab", fem.weights[i], fem.basis, fem.basis)
    np.testing.assert_allclose(full, stock, rtol=1e-11, atol=1e-15)
    min_eigenvalue = np.linalg.eigvalsh(masses).min()
    if min_eigenvalue < -1e-10 * np.abs(masses).max():
        raise RuntimeError("The added drag mass has a nonphysical negative eigenvalue")

    nodes = fem.velocity_dofs[element_ids]
    rows = np.repeat(nodes, 6, axis=1).ravel()
    columns = np.tile(nodes, (1, 6)).ravel()
    drag = alpha * RESISTANCE_SCALE * masses.ravel()
    extra = coo_matrix((np.tile(drag, 2), (np.r_[rows, rows + fem.n_u],
                                          np.r_[columns, columns + fem.n_u])),
                       shape=(fem.size, fem.size)).tocsr()
    return extra, {**diagnostics, "added_penalty_area_mm2": float(extension.area * area_scale),
                   "integrated_added_penalty_area_mm2": integrated_area * area_scale,
                   "elements_with_added_drag": int(len(element_ids)),
                   "subcell_mass_min_eigenvalue": float(min_eigenvalue),
                   "full_element_mass_parity_abs_error": float(np.abs(full - stock).max())}


def run(case, alpha, fitted, *, base_spacing=0.125, target_layer_fraction=0.5,
        max_triangles=170_000, shift_layers=1.):
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run numerical execution as a Sol Slurm job")
    grid = DesignGrid(0.25)
    raw = grid.base_raw if case == "sharp" else obstruction(grid)
    base = build_domain(grid, base_spacing)
    domain, mesh = build_local_binary_domain(base, raw, alpha, target_layer_fraction=target_layer_fraction,
                                             max_triangles=max_triangles)
    extra, correction = assemble_added_drag(domain, grid, raw, alpha, shift_layers)
    resistance = DesignTransfer(domain.cell_index, grid.size).forward(alpha * (1 - raw))

    def solve(fem):
        start = perf_counter()
        state = _continuation(fem, domain.fixed, domain.prescribed, resistance)
        immersed = _baseline_record(domain, state)
        leakage = binary_leakage(domain, state.values, raw)
        split_error = 100 * abs(immersed["side_fraction"] - fitted["side_fraction"])
        power_error = abs(immersed["power_w_per_m"] / fitted["power_w_per_m"] - 1)
        return dict(immersed=immersed, nominal_leakage=leakage, split_error_pp=split_error,
                    power_error_relative=power_error,
                    per_fixture_thresholds_passed=bool(split_error <= 0.02 and power_error <= 0.005
                                                       and leakage["passed"]),
                    elapsed_seconds=perf_counter() - start)

    nominal = solve(domain.fem)
    corrected = copy(domain.fem)
    corrected.viscous_pressure = domain.fem.viscous_pressure + extra
    return dict(case=case, alpha_max=alpha, fitted=fitted, mesh=mesh, correction=correction,
                uncorrected=nominal, corrected=solve(corrected), sol_job_id=os.environ["SLURM_JOB_ID"],
                scope="same mesh, fixed binary engineering diagnostic without production adoption")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("case", choices=("sharp", "obstructed"))
    parser.add_argument("alpha", type=float, choices=PENALTIES)
    parser.add_argument("fitted_reference", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--base-spacing", type=float, default=0.125)
    parser.add_argument("--target-layer-fraction", type=float, default=0.5)
    parser.add_argument("--max-triangles", type=int, default=170_000)
    parser.add_argument("--shift-layers", type=float, default=1.)
    arguments = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical benchmark as a Sol Slurm job")
    with arguments.output.open("x", encoding="utf-8") as stream:
        try:
            fitted = json.loads(arguments.fitted_reference.read_text(encoding="utf-8"))["fixtures"][arguments.case]
            result = run(arguments.case, arguments.alpha, fitted, base_spacing=arguments.base_spacing,
                         target_layer_fraction=arguments.target_layer_fraction,
                         max_triangles=arguments.max_triangles, shift_layers=arguments.shift_layers)
            output = json.dumps({"completed": True, "results": result}, indent=2, allow_nan=False)
        except Exception:
            json.dump({"completed": False, "error": traceback.format_exc()}, stream, indent=2)
            stream.write("\n")
            raise
        stream.write(output + "\n")


if __name__ == "__main__":
    main()
