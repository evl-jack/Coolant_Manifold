"""Simultaneous full-density topology and wall-attached baffle evaluator."""

from dataclasses import replace

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

from ..fem import quadratic_basis
from .density import DensityPipeline, _unit_vector, brinkman, project
from .evaluator import (
    REYNOLDS, RESISTANCE_SCALE, Evaluation, Evaluator, baseline_reference, response_values,
)


BAFFLE_X_MM = (13., 14.)
HEIGHT_MM = 4.
HEIGHT_LIMIT = 3.75 / HEIGHT_MM
TRANSITION_MM = .125


class HybridDensityModel(DensityPipeline):
    """Pin strip inputs to fluid, filter on the full grid, then override projection."""

    def __init__(self, grid, filter_radius_mm=.5):
        super().__init__(grid, filter_radius_mm)
        faces = [(value, grid.x_edges_mm) for value in BAFFLE_X_MM]
        faces.extend((value, grid.y_edges_mm) for value in (0., HEIGHT_MM))
        for value, edges in faces:
            if not np.isclose(edges, value, rtol=0., atol=1e-12).any():
                raise ValueError("The reserved strip must align with design cell faces")
        x, y = grid.centers_mm.T
        self.free = ~((x > BAFFLE_X_MM[0]) & (x < BAFFLE_X_MM[1]) & (y > 0.) & (y < HEIGHT_MM))
        self.size = int(self.free.sum()) + 1
        self.initial = np.r_[grid.base_raw[self.free], 0.]

    def split(self, design):
        design = _unit_vector(design, self.size, "Hybrid design")
        if design[-1] > HEIGHT_LIMIT + 32 * np.finfo(float).eps:
            raise ValueError("Dimensionless baffle height exceeds 3.75/4")
        return design[:-1], float(np.clip(design[-1], 0., HEIGHT_LIMIT))

    def raw(self, design):
        free, _ = self.split(design)
        raw = np.ones(self.grid.size)
        raw[self.free] = free
        return raw

    def evaluate(self, design, q, beta, alpha_max):
        density = super().evaluate(self.raw(design), q, beta, alpha_max)
        projected = density.projected.copy()
        projected[~self.free] = 1.
        projection_derivative = density.projection_derivative.copy()
        projection_derivative[~self.free] = 0.
        resistance, resistance_derivative = brinkman(projected, q, alpha_max)
        area = density.area_mm2 + self.grid.cell_area_mm2 * (projected - density.projected).sum()
        return replace(density, projected=projected, projection_derivative=projection_derivative,
                       resistance=resistance, resistance_derivative=resistance_derivative,
                       area_mm2=float(area))

    def project(self, design, beta=0.):
        filtered = self.filter(self.raw(design))
        projected, _ = project(filtered, beta, self.threshold)
        projected[~self.free] = 1.
        return projected

    def projected_pullback(self, density, covector):
        return self.filter_pullback(density.projection_derivative * covector)[self.free]


def _gauss(order):
    points, weights = np.polynomial.legendre.leggauss(order)
    return (points + 1.) / 2., weights / 2.


class BaffleSubcells:
    """Integrate horizontal pieces of each parent P2 element in the strip."""

    def __init__(self, domain, model, transition_mm=TRANSITION_MM):
        if not np.isfinite(transition_mm) or transition_mm <= 0:
            raise ValueError("Baffle transition length must be finite and positive")
        self.transition_mm = float(transition_mm)
        self.fem = domain.fem
        design = np.flatnonzero(domain.cell_index >= 0)
        element_ids = design[~model.free[domain.cell_index[design]]]
        self.elements = element_ids
        vertices = self.fem.points[self.fem.triangles[element_ids]] * HEIGHT_MM
        self.anchor = vertices[:, 0]
        self.inverse = np.linalg.inv((vertices[:, 1:] - vertices[:, :1]).transpose(0, 2, 1))
        levels = np.sort(vertices[..., 1], axis=1)
        parent, layer = np.nonzero(np.diff(levels, axis=1) > 0.)
        self.parent = parent
        self.lower, self.upper = levels[parent, layer], levels[parent, layer + 1]
        a, b = vertices[parent], np.roll(vertices[parent], -1, axis=1)
        middle = (self.lower + self.upper) / 2.
        crosses = ((np.minimum(a[..., 1], b[..., 1]) < middle[:, None])
                   & (middle[:, None] < np.maximum(a[..., 1], b[..., 1])))
        a, b = a[crosses].reshape(-1, 2, 2), b[crosses].reshape(-1, 2, 2)
        slope = (b[..., 0] - a[..., 0]) / (b[..., 1] - a[..., 1])
        intercept = a[..., 0] - slope * a[..., 1]
        order = np.argsort(intercept + slope * middle[:, None], axis=1)
        self.slope = np.take_along_axis(slope, order, axis=1)
        self.intercept = np.take_along_axis(intercept, order, axis=1)
        self.y_nodes, self.y_weights = _gauss(6)
        self.x_nodes, self.x_weights = _gauss(3)

        nodes = self.fem.velocity_dofs[element_ids]
        rows = np.repeat(nodes, 6, axis=1).ravel()
        columns = np.tile(nodes, (1, 6)).ravel()
        self.rows = np.r_[rows, rows + self.fem.n_u]
        self.columns = np.r_[columns, columns + self.fem.n_u]

    def _integrate(self, lower, upper, height, transition):
        active = np.flatnonzero(upper > lower)
        parent = self.parent[active]
        if not active.size:
            empty = np.empty((0, 6, 6))
            return parent, empty, empty
        span = upper[active] - lower[active]
        y = lower[active, None] + span[:, None] * self.y_nodes
        limits = self.intercept[active, None] + self.slope[active, None] * y[..., None]
        widths = limits[..., 1] - limits[..., 0]
        x = limits[..., :1] + widths[..., None] * self.x_nodes
        points = np.stack((x, np.broadcast_to(y[..., None], x.shape)), axis=-1)
        xy = np.einsum("eyxj,eij->eyxi", points - self.anchor[parent, None, None],
                       self.inverse[parent])
        bary = np.concatenate((1. - xy.sum(axis=-1, keepdims=True), xy), axis=-1)
        basis, _ = quadratic_basis(bary)
        weights = (span[:, None, None] * widths[..., None] * self.y_weights[None, :, None]
                   * self.x_weights / HEIGHT_MM**2)
        if transition:
            t = (height - y) / self.transition_mm
            gate = t**3 * (10. + t * (-15. + 6. * t))
            derivative = (30. * HEIGHT_MM / self.transition_mm) * t**2 * (1. - t)**2
        else:
            gate, derivative = np.ones_like(y), np.zeros_like(y)
        mass = np.einsum("eyx,ey,eyxi,eyxj->eij", weights, gate, basis, basis, optimize=True)
        tangent = np.einsum("eyx,ey,eyxi,eyxj->eij", weights, derivative, basis, basis, optimize=True)
        return parent, mass, tangent

    def integrate(self, s, alpha_max):
        height = HEIGHT_MM * s
        shoulder = height - self.transition_mm
        mass = np.zeros((len(self.elements), 6, 6))
        tangent = np.zeros_like(mass)
        for lower, upper, transition in (
            (self.lower, np.minimum(self.upper, shoulder), False),
            (np.maximum(self.lower, shoulder), np.minimum(self.upper, height), True),
        ):
            parent, addition, derivative = self._integrate(lower, upper, height, transition)
            np.add.at(mass, parent, addition)
            np.add.at(tangent, parent, derivative)

        def matrix(local):
            entries = RESISTANCE_SCALE * alpha_max * local.ravel()
            return coo_matrix((np.tile(entries, 2), (self.rows, self.columns)),
                              shape=(self.fem.size, self.fem.size)).tocsr()

        scale = HEIGHT_MM**2
        return matrix(mass), matrix(tangent), scale * mass.sum(), scale * tangent.sum()


class HybridEvaluator(Evaluator):
    """Reuse the Navier–Stokes state and adjoint with integrated baffle drag."""

    def __init__(self, domain, model, stage, power_reference, *, initial_state=None,
                 height_cost=0., graph=None):
        super().__init__(domain, model, stage, power_reference, initial_state=initial_state)
        domain.fem.inertia = REYNOLDS
        self.baffle = BaffleSubcells(domain, model)
        self.graph = graph
        self.height_cost = height_cost
        self._operator_design = None
        self._operator_data = None

    def operators(self, design):
        if self._operator_design is not None and np.array_equal(design, self._operator_design):
            return self._operator_data
        density = self.density(design)
        _, s = self.pipeline.split(design)
        extra, tangent, solid_area, solid_gradient = self.baffle.integrate(s, self.stage.alpha_max)
        drag = self.domain.fem.drag_matrix(self.transfer.forward(density.resistance)) + extra
        area_gradient = np.r_[self.pipeline.area_gradient(density)[self.pipeline.free], -solid_gradient]
        self._operator_design = np.array(design, copy=True)
        self._operator_data = (drag, tangent, density.area_mm2 - solid_area, area_gradient)
        return self._operator_data

    def residual_direction(self, design, state_values, direction):
        density = self.density(design)
        raw_direction = np.zeros(self.pipeline.grid.size)
        raw_direction[self.pipeline.free] = direction[:-1]
        resistance_direction = (density.resistance_derivative * density.projection_derivative
                                * self.pipeline.filter(raw_direction))
        topology = self.domain.fem.drag_matrix(self.transfer.forward(resistance_direction))
        _, tangent, _, _ = self.operators(design)
        return (topology + direction[-1] * tangent) @ state_values

    def certificate(self, design):
        if self.graph is None:
            raise ValueError("No connectivity graph was configured")
        return self.graph.evaluate(self.density(design).projected, baffle_height_mm=HEIGHT_MM * design[-1])

    def evaluate(self, design):
        design = np.asarray(design, dtype=float)
        if self.last_design is not None and np.array_equal(design, self.last_design):
            return self.last_evaluation
        density = self.density(design)
        drag, height_tangent, area, area_gradient = self.operators(design)
        fem = self.domain.fem
        self.forward_attempts += 1
        state = fem.solve(self.domain.fixed, self.domain.prescribed, self.load, added_drag=drag,
                          initial=None if self.state is None else self.state.values)
        self.forward_solves += 1
        self.newton_iterations += len(state.residual_history) - 1
        flux, imbalance, power, fraction, state_gradients = response_values(self.domain, state.values)
        free = state.free
        tangent = splu(state.tangent[free][:, free].tocsc())
        adjoints = np.zeros((fem.size, 2))
        adjoints[free] = tangent.solve(state_gradients[:, free].T, trans="T")
        self.adjoint_rhs += 2
        velocity = fem.velocity_fields(state.values)[0]
        height_action = height_tangent @ state.values
        gradients = []
        for adjoint in adjoints.T:
            adjoint_velocity = fem.velocity_fields(adjoint)[0]
            covector = -fem.weights * np.sum(velocity * adjoint_velocity, axis=2)
            topology = self.pipeline.pullback(density, self.transfer.transpose(covector))
            gradients.append(np.r_[topology[self.pipeline.free], -adjoint @ height_action])
        power_gradient, fraction_gradient = gradients
        objective_gradient = power_gradient / self.power_reference
        objective_gradient[-1] += self.height_cost
        constraints = np.array([(0.499 - fraction) / .001, (fraction - .501) / .001, area - 83.])
        jacobian = np.stack((-fraction_gradient / .001, fraction_gradient / .001, area_gradient))
        if self.graph is not None:
            certificate = self.certificate(design)
            scale = self.graph.minimum_conductance
            graph_topology = np.stack([self.pipeline.projected_pullback(density, row)
                                       for row in certificate.constraint_projected_jacobian])
            graph_height = HEIGHT_MM * certificate.constraint_height_jacobian
            constraints = np.r_[constraints, certificate.constraints / scale]
            jacobian = np.vstack((jacobian, np.column_stack((graph_topology, graph_height)) / scale))
        result = Evaluation(power / self.power_reference + self.height_cost * design[-1],
                            objective_gradient, constraints, jacobian, power, fraction, area, flux,
                            imbalance, state.residual_history[-1],
                            float(4. * np.mean(density.projected * (1. - density.projected))))
        self.state, self.last_design, self.last_evaluation = state, design.copy(), result
        return result


def hybrid_reference(domain, model, alpha_max):
    """Legacy sharp-seeded zero-height reference, with the original filter."""
    pipeline = DensityPipeline(model.grid, model.filter_radius_mm, model.threshold)
    return baseline_reference(domain, pipeline, alpha_max)
