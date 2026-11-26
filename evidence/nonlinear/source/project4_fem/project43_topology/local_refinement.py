"""Conforming local analysis meshes for the fixed, unfiltered Gate 3 fixtures."""

import numpy as np

from ..fem import LOCAL_EDGES, TaylorHood
from ..junction import boundary_data
from .domain import REFERENCE_HEIGHT_MM, TopologyDomain
from .evaluator import VISCOSITY


def _edge_topology(triangles):
    pairs = np.sort(triangles[:, LOCAL_EDGES], axis=-1).reshape(-1, 2)
    edges, inverse, incidence = np.unique(pairs, axis=0, return_inverse=True, return_counts=True)
    return edges, inverse.reshape(-1, 3), incidence


def _merged_axes(segments):
    """Merge collinear horizontal or vertical segments into closed intervals."""
    runs = {}
    for start, end in segments:
        axes = np.flatnonzero(start != end)
        if len(axes) != 1:
            raise ValueError("Binary interface and port edges must be axis aligned")
        axis = int(axes[0])
        level = start[1-axis]
        runs.setdefault((axis, level), []).append(sorted((start[axis], end[axis])))

    merged = []
    for (axis, level), intervals in sorted(runs.items()):
        ordered = sorted(intervals)
        first, last = ordered[0]
        for low, high in ordered[1:]:
            if low <= last:
                last = max(last, high)
            else:
                merged.append((axis, level, first, last))
                first, last = low, high
        merged.append((axis, level, first, last))
    return merged


def _measure(vertices, segments, band_mm):
    """Use overlapping bounding boxes to select a superset of band triangles."""
    lower, upper = vertices.min(axis=1), vertices.max(axis=1)
    selected = np.zeros(len(vertices), dtype=bool)
    normal_spacing, tangential_spacing = 0.0, 0.0
    for axis, level, first, last in segments:
        normal = 1-axis
        in_band = ((lower[:, normal] <= level+band_mm) & (upper[:, normal] >= level-band_mm)
                   & (lower[:, axis] <= last+band_mm) & (upper[:, axis] >= first-band_mm))
        selected |= in_band
        if np.any(in_band):
            span = upper[in_band] - lower[in_band]
            normal_spacing = max(normal_spacing, float(span[:, normal].max()))
            tangential_spacing = max(tangential_spacing, float(span[:, axis].max()))
    return selected, normal_spacing, tangential_spacing


def _split_edges(points, triangles, parents, edges, element_edges, marked, maximum_triangles):
    """Use the same midpoint on both sides of every split edge."""
    midpoints = np.full(len(edges), -1, dtype=np.intp)
    midpoints[marked] = len(points) + np.arange(np.count_nonzero(marked))
    element_mids = midpoints[element_edges]
    counts = np.count_nonzero(element_mids >= 0, axis=1)
    total = len(triangles) + counts.sum()
    if total > maximum_triangles:
        raise RuntimeError(f"Local refinement would need {total} triangles, above the cap {maximum_triangles}")
    points = np.vstack((points, points[edges[marked]].mean(axis=1)))
    refined, source = [triangles[counts == 0]], [parents[counts == 0]]
    for count in (1, 2, 3):
        keep = counts == count
        if not np.any(keep):
            continue
        vertices, mids = triangles[keep], element_mids[keep]
        if count != 3:
            # For one split, rotate its edge into slot 0. For two splits,
            # rotate the unsplit edge into slot 2. Cyclic rotations preserve orientation.
            pivot = np.argmax(mids >= 0, axis=1) if count == 1 else (np.argmin(mids >= 0, axis=1) + 1) % 3
            order = (pivot[:, None] + np.arange(3)) % 3
            vertices = np.take_along_axis(vertices, order, axis=1)
            mids = np.take_along_axis(mids, order, axis=1)
        a, b, c = vertices.T
        m, n, p = mids.T
        if count == 1:
            children = ((a, m, c), (m, b, c))
        elif count == 2:
            # The corner at b leaves a quadrilateral. Use its shorter diagonal.
            ac = np.linalg.norm(points[a] - points[n], axis=1)
            cm = np.linalg.norm(points[c] - points[m], axis=1)
            use_an = ac <= cm
            children = ((m, b, n), (a, m, np.where(use_an, n, c)),
                        (np.where(use_an, a, m), n, c))
        else:
            children = ((a, m, p), (m, b, n), (p, n, c), (m, n, p))
        refined += [np.column_stack(child) for child in children]
        source += [parents[keep]] * len(children)
    return points, np.vstack(refined), np.concatenate(source)


def _double_area(vertices):
    ab, ac = vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]
    return ab[:, 0]*ac[:, 1] - ab[:, 1]*ac[:, 0]


def build_local_binary_domain(domain, raw, alpha_max, *, target_layer_fraction=0.5,
                              band_layer_multiple=2.0, port_band_mm=0.125,
                              port_max_diameter_mm=0.0625, max_triangles=240_000, max_passes=24):
    """Refine a single fixed binary fixture and return its domain plus diagnostics.

    Mark a conservatively covered rectangular band on both sides of each
    original solid/fluid segment. Bound the triangle's normal and tangential
    spans, including near endpoints and corners of the axis-aligned segments.
    A separate, wider-spacing band covers all through-flow ports and their
    endpoints. Fluid walls away from those bands retain the passively refined
    base mesh. This certificate covers the supplied unfiltered binary fixture.
    Moving or projected interfaces require a separate analysis mesh assessment.
    """
    if (not np.isfinite([alpha_max, target_layer_fraction, band_layer_multiple,
                         port_band_mm, port_max_diameter_mm]).all()
            or min(alpha_max, target_layer_fraction, band_layer_multiple,
                   port_band_mm, port_max_diameter_mm) <= 0):
        raise ValueError("Resistance and refinement lengths must be finite and positive")
    if max_triangles < len(domain.fem.triangles) or max_passes < 1:
        raise ValueError("Refinement limits must permit the base mesh and at least one pass")
    raw = np.asarray(raw)
    indices = np.asarray(domain.cell_index)
    if raw.ndim != 1 or not np.all((raw == 0) | (raw == 1)) or np.any(indices >= len(raw)):
        raise ValueError("Provide one binary raw value per design cell")
    designable = indices >= 0
    solid = np.zeros(len(indices), dtype=bool)
    solid[designable] = raw[indices[designable]] == 0

    points = domain.fem.points * REFERENCE_HEIGHT_MM
    triangles = domain.fem.triangles
    original_area = _double_area(points[triangles])
    edges, element_edges, incidence = _edge_topology(triangles)
    solid_incidence = np.bincount(element_edges[solid].ravel(), minlength=len(edges))
    interface = (incidence == 2) & (solid_incidence == 1)
    if not np.any(interface):
        raise ValueError("The fixture has no interior binary fluid/solid interface")
    segments = _merged_axes(points[edges[interface]])
    port_segments = _merged_axes(points[np.vstack(domain.ports)[:, :2]])

    layer_mm = 1000.0*np.sqrt(VISCOSITY/alpha_max)
    band_mm = band_layer_multiple*layer_mm
    target_mm = target_layer_fraction*layer_mm
    parents = np.arange(len(triangles), dtype=np.intp)
    counts = [len(triangles)]
    for pass_number in range(max_passes + 1):
        vertices = points[triangles]
        diameters = np.linalg.norm(vertices[:, [1, 2, 0]] - vertices, axis=2).max(axis=1)
        axis_span = np.ptp(vertices, axis=1).max(axis=1)
        in_interface_band, normal_mm, tangent_mm = _measure(vertices, segments, band_mm)
        in_port_band, port_normal_mm, port_tangent_mm = _measure(vertices, port_segments, port_band_mm)
        must_refine = ((in_interface_band & (axis_span > target_mm))
                       | (in_port_band & (diameters > port_max_diameter_mm)))
        if not np.any(must_refine):
            break
        if pass_number == max_passes:
            raise RuntimeError("Local refinement exhausted its pass budget before meeting the spacing bounds")
        marked = np.zeros(len(edges), dtype=bool)
        marked[element_edges[must_refine].ravel()] = True
        points, triangles, parents = _split_edges(points, triangles, parents, edges, element_edges,
                                                   marked, max_triangles)
        counts.append(len(triangles))
        edges, element_edges, incidence = _edge_topology(triangles)

    refined_area = _double_area(vertices)
    max_squared_edge = diameters**2
    minimum_area_quality = float(np.min(refined_area / max_squared_edge))
    if minimum_area_quality < 1e-3 or np.any(incidence > 2):
        raise RuntimeError("Local refinement produced a degenerate or nonmanifold mesh")
    np.testing.assert_allclose(np.bincount(parents, weights=refined_area, minlength=len(original_area)),
                               original_area, rtol=1e-12, atol=1e-14)
    old_boundary = domain.fem.boundary_edges[:, :2]
    old_vertices = domain.fem.points[old_boundary] * REFERENCE_HEIGHT_MM
    boundary_vertices = points[edges[incidence == 1]]
    # A hanging edge inside a flat mesh would add to topological boundary length.
    old_length = np.linalg.norm(old_vertices[:, 1] - old_vertices[:, 0], axis=1).sum()
    final_length = np.linalg.norm(boundary_vertices[:, 1] - boundary_vertices[:, 0], axis=1).sum()
    np.testing.assert_allclose(final_length, old_length, rtol=1e-12)

    fem = TaylorHood(points / REFERENCE_HEIGHT_MM, triangles, viscosity=domain.fem.viscosity,
                     inertia=domain.fem.inertia)
    fixed, prescribed, ports, flux = boundary_data(fem)
    for old_port, final_port in zip(domain.ports, ports):
        old_vertices = domain.fem.points[old_port[:, :2]]
        new_vertices = fem.points[final_port[:, :2]]
        np.testing.assert_allclose(new_vertices.min(axis=(0, 1)), old_vertices.min(axis=(0, 1)), atol=1e-14)
        np.testing.assert_allclose(new_vertices.max(axis=(0, 1)), old_vertices.max(axis=(0, 1)), atol=1e-14)
        np.testing.assert_allclose(np.linalg.norm(new_vertices[:, 1] - new_vertices[:, 0], axis=1).sum(),
                                   np.linalg.norm(old_vertices[:, 1] - old_vertices[:, 0], axis=1).sum(), rtol=1e-12)
    refined = TopologyDomain(fem, fixed, prescribed, ports, flux, indices[parents], domain.max_spacing_mm)
    diagnostics = {
        "scope": "fixed unfiltered binary fixture only",
        "brinkman_length_mm": layer_mm,
        "interface_band_half_width_mm": band_mm,
        "interface_normal_and_tangential_target_mm": target_mm,
        "interface_max_diameter_mm": float(diameters[in_interface_band].max()),
        "interface_max_normal_mm": normal_mm,
        "interface_max_tangential_mm": tangent_mm,
        "interface_segment_count": len(segments),
        "port_band_half_width_mm": port_band_mm,
        "port_diameter_target_mm": port_max_diameter_mm,
        "port_max_diameter_mm": float(diameters[in_port_band].max()),
        "port_max_normal_mm": port_normal_mm,
        "port_max_tangential_mm": port_tangent_mm,
        "minimum_double_area_to_diameter_squared": minimum_area_quality,
        "triangles_by_pass": counts,
        "dofs": fem.size,
    }
    return refined, diagnostics
