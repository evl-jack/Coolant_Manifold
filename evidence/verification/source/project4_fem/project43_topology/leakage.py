"""P2 leakage diagnostics for Gate 3 element-aligned binary solid regions."""

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from ..fem import LOCAL_EDGES
from ..project43_leakage import directional_integrals
from .domain import REFERENCE_HEIGHT_MM


GROSS_LIMIT = .001
SECTION_ABSOLUTE_LIMIT = .0002
SOLID_RMS_LIMIT = .01
_SECTIONS = (
    ("roof_x8", 0, 8., 4., 6.),
    ("roof_x12", 0, 12., 4., 6.),
    ("roof_x14", 0, 14., 4., 6.),
    ("side_shortcut_y5", 1, 5., 6., 16.),
)


def _edge_topology(fem):
    """Reuse each P2 midpoint DOF as the unique index of its mesh edge."""
    count = fem.n_u - fem.n_p
    first = np.full(count, len(fem.triangles), dtype=np.intp)
    last = np.full(count, -1, dtype=np.intp)
    vertices = np.empty((count, 2), dtype=np.intp)
    owner = np.arange(len(fem.triangles))
    for slot, (i, j) in enumerate(LOCAL_EDGES):
        edge = fem.velocity_dofs[:, 3 + slot] - fem.n_p
        np.minimum.at(first, edge, owner)
        np.maximum.at(last, edge, owner)
        vertices[edge] = np.sort(fem.triangles[:, [i, j]], axis=1)
    return vertices, first, last


def _edge_flux(edges, vertices, points, velocity, n_p, normals):
    start, end = vertices[edges].T
    lengths = np.linalg.norm(points[end] - points[start], axis=1)
    nodes = np.column_stack((start, n_p + edges, end))
    weighted_normal = lengths[:, None] * np.broadcast_to(normals, (len(edges), 2))
    values = np.einsum("evi,ei->ev", velocity[nodes], weighted_normal)
    totals = np.empty((len(edges), 2))
    for i, sample in enumerate(values):
        totals[i] = directional_integrals(*sample)
    return lengths, totals


def _solid_statistics(fem, velocity, solid_ids, component, count):
    component_area = np.zeros(count)
    sum_speed2 = np.zeros(count)
    peak = np.zeros(count)
    for start in range(0, len(solid_ids), 4096):
        block = slice(start, start + 4096)
        elements = solid_ids[block]
        labels = component[block]
        nodal = velocity[fem.velocity_dofs[elements]]
        quadrature_velocity = np.einsum("qa,eai->eqi", fem.basis, nodal)
        weights = fem.weights[elements]
        speed2 = np.sum(quadrature_velocity ** 2, axis=2)
        np.add.at(component_area, labels, weights.sum(axis=1))
        np.add.at(sum_speed2, labels, np.sum(weights * speed2, axis=1))
        maximum = np.maximum(speed2.max(axis=1), np.sum(nodal ** 2, axis=2).max(axis=1))
        np.maximum.at(peak, labels, maximum)
    return component_area, sum_speed2, np.sqrt(peak)


def _section_report(vertices, fem, velocity, both_solid, inlet):
    def matches(values, expected):
        return np.isclose(values, expected, rtol=0., atol=1e-11)

    endpoints = fem.points[vertices] * REFERENCE_HEIGHT_MM
    reports = {}
    for name, axis, location, lower, upper in _SECTIONS:
        on_line = np.all(matches(endpoints[:, :, axis], location), axis=1)
        tangent_endpoints = endpoints[:, :, 1 - axis]
        inside = np.all((tangent_endpoints >= lower - 1e-11) & (tangent_endpoints <= upper + 1e-11), axis=1)
        represented = on_line & inside
        represented_length = np.abs(np.diff(tangent_endpoints[represented], axis=1)).sum()
        if not matches(represented_length, upper - lower):
            raise ValueError(f"The FEM mesh does not represent the complete {name} leakage cut")
        edges = np.flatnonzero(represented & both_solid)
        normal = np.eye(2)[axis]
        lengths, totals = _edge_flux(edges, vertices, fem.points, velocity, fem.n_p, normal)
        positive, negative = totals.sum(axis=0) / inlet
        absolute = positive + negative
        reports[name] = dict(normal=normal.tolist(), cut_length_mm=float(represented_length),
                             solid_length_mm=float(lengths.sum() * REFERENCE_HEIGHT_MM),
                             positive_fraction=float(positive), negative_fraction=float(negative),
                             signed_fraction=float(positive - negative), absolute_fraction=float(absolute),
                             passed=bool(absolute <= SECTION_ABSOLUTE_LIMIT))
    return reports


def binary_leakage(domain, state_values, raw_density):
    """Measure gross interface crossings and internal solid transport.

    Velocity is in inlet bulk-speed units. All fluxes are fractions of the
    signed inlet flow. A through-solid cut includes traces with solid on both
    sides, excluding the fluid-solid interfaces themselves.
    """
    fem = domain.fem
    raw_density = np.asarray(raw_density, dtype=float)
    state_values = np.asarray(state_values, dtype=float)
    if state_values.shape != (fem.size,) or not np.isfinite(state_values).all():
        raise ValueError("Leakage calculation requires a finite full FEM state vector")
    inside = domain.cell_index >= 0
    if (raw_density.shape != (int(domain.cell_index.max()) + 1,)
            or not np.all((raw_density == 0.) | (raw_density == 1.))):
        raise ValueError("Binary leakage calculation requires a 0/1 density for each design cell")
    inlet = float(domain.flux[0] @ state_values)
    if inlet <= 0.:
        raise ValueError("Inlet flux must be positive for normalized leakage calculations")
    velocity = state_values[:2 * fem.n_u].reshape(2, -1).T
    solid = np.zeros(len(fem.triangles), dtype=bool)
    solid[inside] = raw_density[domain.cell_index[inside]] == 0.
    solid_ids = np.flatnonzero(solid)
    if not len(solid_ids):
        raise ValueError("Gate 3 leakage fixtures must contain solid material")

    vertices, first, last = _edge_topology(fem)
    internal = first != last
    both_solid = internal & solid[first] & solid[last]
    compact = np.cumsum(solid) - 1
    graph_i = compact[first[both_solid]]
    graph_j = compact[last[both_solid]]
    adjacency = coo_matrix((np.ones(len(graph_i)), (graph_i, graph_j)), shape=(len(solid_ids), len(solid_ids)))
    count, component = connected_components(adjacency, directed=False)
    area, speed2, peak = _solid_statistics(fem, velocity, solid_ids, component, count)

    interface = np.flatnonzero(internal & (solid[first] != solid[last]))
    solid_owner = np.where(solid[first[interface]], first[interface], last[interface])
    a, b = fem.points[vertices[interface]].transpose(1, 0, 2)
    delta = b - a
    normals = np.column_stack((delta[:, 1], -delta[:, 0]))
    normals /= np.linalg.norm(delta, axis=1)[:, None]
    solid_center = fem.points[fem.triangles[solid_owner]].mean(axis=1)
    inward = np.sum(normals * (solid_center - a), axis=1) > 0.
    normals[~inward] *= -1.
    lengths, exchange = _edge_flux(interface, vertices, fem.points, velocity, fem.n_p, normals)
    labels = component[compact[solid_owner]]
    component_lengths = np.bincount(labels, weights=lengths, minlength=count) * REFERENCE_HEIGHT_MM
    gross_by_component = np.zeros((count, 2))
    np.add.at(gross_by_component, labels, exchange / inlet)
    gross_ingress, gross_egress = gross_by_component.sum(axis=0)
    rms = float(np.sqrt(speed2.sum() / area.sum()))
    sections = _section_report(vertices, fem, velocity, both_solid, inlet)
    components = []
    for i, (ingress, egress) in enumerate(gross_by_component):
        components.append(dict(component=i, area_mm2=float(area[i] * REFERENCE_HEIGHT_MM ** 2),
                               interface_length_mm=float(component_lengths[i]),
                               ingress_fraction=float(ingress), egress_fraction=float(egress),
                               net_ingress_fraction=float(ingress - egress),
                               rms_speed_over_bulk=float(np.sqrt(speed2[i] / area[i])),
                               sampled_max_speed_over_bulk=float(peak[i])))

    checks = dict(gross_ingress=bool(gross_ingress <= GROSS_LIMIT),
                  gross_egress=bool(gross_egress <= GROSS_LIMIT),
                  sections=all(report["passed"] for report in sections.values()),
                  solid_rms=bool(rms <= SOLID_RMS_LIMIT))
    return dict(passed=all(checks.values()), checks=checks, inlet_flow_dimensionless=inlet,
                limits=dict(gross_each=GROSS_LIMIT, section_absolute_each=SECTION_ABSOLUTE_LIMIT,
                            solid_rms_over_bulk=SOLID_RMS_LIMIT),
                gross_ingress_fraction=float(gross_ingress), gross_egress_fraction=float(gross_egress),
                net_ingress_fraction=float(gross_ingress - gross_egress),
                interface_length_mm=float(component_lengths.sum()), area_mm2=float(area.sum() * REFERENCE_HEIGHT_MM ** 2),
                rms_speed_over_bulk=rms, sampled_max_speed_over_bulk=float(peak.max()),
                sections=sections, solid_components=components)
