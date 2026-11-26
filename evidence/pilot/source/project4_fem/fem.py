"""Self-written triangular P2/P1 finite elements for incompressible flow."""

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


BARYCENTRIC_GRADIENTS = np.array([[-1., -1.], [1., 0.], [0., 1.]])
LOCAL_EDGES = ((0, 1), (1, 2), (2, 0))


def triangle_rule(order=5):
    """Tensor Gauss rule mapped onto the unit triangle, whose area is 1/2."""
    points, weights = np.polynomial.legendre.leggauss(order)
    points, weights = (points + 1) / 2, weights / 2
    r, s = np.meshgrid(points, points, indexing="ij")
    x, y = r.ravel(), ((1 - r) * s).ravel()
    return np.column_stack((1 - x - y, x, y)), (np.outer(weights, weights) * (1 - r)).ravel()


def quadratic_basis(barycentric):
    """Vertex functions followed by midpoint functions on edges 01, 12, 20."""
    values = [barycentric * (2 * barycentric - 1)]
    gradients = [(4 * barycentric - 1)[..., None] * BARYCENTRIC_GRADIENTS]
    for i, j in LOCAL_EDGES:
        a, b = barycentric[..., i], barycentric[..., j]
        values.append((4 * a * b)[..., None])
        gradients.append((4 * (a[..., None] * BARYCENTRIC_GRADIENTS[j]
                               + b[..., None] * BARYCENTRIC_GRADIENTS[i]))[..., None, :])
    return np.concatenate(values, axis=-1), np.concatenate(gradients, axis=-2)


def rectangle(nx, ny, length=1., height=1.):
    x, y = np.meshgrid(np.linspace(0., length, nx + 1), np.linspace(0., height, ny + 1), indexing="xy")
    points = np.column_stack((x.ravel(), y.ravel()))
    nodes = np.arange(len(points)).reshape(ny + 1, nx + 1)
    a, b, c, d = [corner.ravel() for corner in
                  (nodes[:-1, :-1], nodes[:-1, 1:], nodes[1:, 1:], nodes[1:, :-1])]
    return points, np.vstack((np.column_stack((a, b, c)), np.column_stack((a, c, d))))


@dataclass
class State:
    values: np.ndarray
    residual_history: list
    tangent: object
    free: np.ndarray


class TaylorHood:
    def __init__(self, points, triangles, *, viscosity=1., inertia=0.):
        self.points = np.asarray(points, dtype=float)
        self.triangles = np.asarray(triangles, dtype=int)
        self.viscosity, self.inertia = viscosity, inertia
        vertices = self.points[self.triangles]
        jacobian = (vertices[:, 1:] - vertices[:, :1]).transpose(0, 2, 1)
        determinant = np.linalg.det(jacobian)
        if np.any(determinant <= 0):
            raise ValueError("Triangles must have positive orientation and area")

        edge_ends = np.sort(self.triangles[:, LOCAL_EDGES], axis=-1).reshape(-1, 2)
        edges, incidence, count = np.unique(edge_ends, axis=0, return_inverse=True, return_counts=True)
        self.n_p = len(self.points)
        self.velocity_points = np.vstack((self.points, self.points[edges].mean(axis=1)))
        self.n_u = len(self.velocity_points)
        self.velocity_dofs = np.column_stack((self.triangles, self.n_p + incidence.reshape(-1, 3)))
        self.boundary_edges = np.c_[edges[count == 1], self.n_p + np.flatnonzero(count == 1)]
        self.dofs = np.column_stack((self.velocity_dofs, self.n_u + self.velocity_dofs,
                                     2 * self.n_u + self.triangles))
        self.size = 2 * self.n_u + self.n_p
        self.rows = np.broadcast_to(self.dofs[:, :, None], (len(vertices), 15, 15)).ravel()
        self.cols = np.broadcast_to(self.dofs[:, None, :], (len(vertices), 15, 15)).ravel()

        self.p_basis, weights = triangle_rule()
        self.basis, reference_gradient = quadratic_basis(self.p_basis)
        self.gradient = np.einsum("qai,eij->eqaj", reference_gradient, np.linalg.inv(jacobian))
        self.weights = determinant[:, None] * weights
        self.quad_points = np.einsum("qa,eai->eqi", self.p_basis, vertices)
        self.viscous_pressure = self.assemble_linear()

    def sparse(self, local):
        return coo_matrix((local.ravel(), (self.rows, self.cols)), shape=(self.size, self.size)).tocsr()

    def assemble_linear(self):
        local = np.zeros((len(self.triangles), 15, 15))
        w, g, p = self.weights, self.gradient, self.p_basis
        grad_dot = np.einsum("eq,eqai,eqbi->eab", w, g, g)
        for i in range(2):
            a = slice(6 * i, 6 * (i + 1))
            coupling = -np.einsum("eq,qa,eqb->eab", w, p, g[..., i])
            local[:, 12:, a], local[:, a, 12:] = coupling, coupling.transpose(0, 2, 1)
            for j in range(2):
                b = slice(6 * j, 6 * (j + 1))
                cross = np.einsum("eq,eqa,eqb->eab", w, g[..., j], g[..., i])
                local[:, a, b] = self.viscosity * (cross + (grad_dot if i == j else 0.))
        return self.sparse(local)

    def drag_matrix(self, resistance):
        w = self.weights * np.asarray(resistance)
        mass = np.einsum("eq,qa,qb->eab", w, self.basis, self.basis)
        local = np.zeros((len(self.triangles), 15, 15))
        local[:, :6, :6] = local[:, 6:12, 6:12] = mass
        return self.sparse(local)

    def velocity_fields(self, state):
        nodal = np.c_[state[:self.n_u], state[self.n_u:2 * self.n_u]][self.velocity_dofs]
        return (np.einsum("qa,eai->eqi", self.basis, nodal),
                np.einsum("eqaj,eai->eqij", self.gradient, nodal))

    def scatter_velocity(self, values):
        rhs = np.zeros(self.size)
        for i in range(2):
            np.add.at(rhs, self.velocity_dofs + i * self.n_u, values[..., i])
        return rhs

    def load(self, force=None, traction=None):
        values = np.zeros(self.size)
        if force is not None:
            body = force(self.quad_points)
            values = self.scatter_velocity(np.einsum("eq,qa,eqi->eai", self.weights, self.basis, body))
        if traction is not None:
            edges, function = traction
            points, weights = np.polynomial.legendre.leggauss(4)
            s, w = (points + 1) / 2, weights / 2
            start, stop = self.points[edges[:, 0]], self.points[edges[:, 1]]
            sample = start[:, None] + s[None, :, None] * (stop - start)[:, None]
            shape = np.column_stack(((1 - s) * (1 - 2 * s), s * (2 * s - 1), 4 * s * (1 - s)))
            line = np.einsum("e,q,qa,eqi->eai", np.linalg.norm(stop - start, axis=1), w, shape, function(sample))
            for i in range(2):
                np.add.at(values, edges + i * self.n_u, line[..., i])
        return values

    def residual(self, state, operator, load, tangent=True):
        value = operator @ state - load
        if not self.inertia:
            return value, operator if tangent else None
        velocity, jac_velocity = self.velocity_fields(state)
        advection = np.einsum("eqj,eqij->eqi", velocity, jac_velocity)
        w, n = self.inertia * self.weights, self.basis
        value += self.scatter_velocity(np.einsum("eq,qa,eqi->eai", w, n, advection))
        if not tangent:
            return value, None
        directional = np.einsum("eqj,eqbj->eqb", velocity, self.gradient)
        convected = np.einsum("eq,qa,eqb->eab", w, n, directional)
        local = np.zeros((len(self.triangles), 15, 15))
        for i in range(2):
            for j in range(2):
                transported = np.einsum("eq,qa,qb,eq->eab", w, n, n, jac_velocity[..., i, j])
                local[:, 6*i:6*(i+1), 6*j:6*(j+1)] = transported + (convected if i == j else 0.)
        return value, operator + self.sparse(local)

    def solve(self, fixed, prescribed, load, *, resistance=0., initial=None, tolerance=1e-11,
              iterations=35, added_drag=None):
        operator = self.viscous_pressure + self.drag_matrix(resistance)
        if added_drag is not None:
            operator = operator + added_drag
        free = np.setdiff1d(np.arange(self.size), fixed)
        state = np.zeros(self.size) if initial is None else initial.copy()
        state[fixed] = prescribed
        history = []
        for _ in range(iterations + 1):
            residual, jacobian = self.residual(state, operator, load)
            error = np.linalg.norm(residual[free], ord=np.inf)
            history.append(float(error))
            if error <= tolerance:
                return State(state, history, jacobian, free)
            update = splu(jacobian[free][:, free].tocsc()).solve(-residual[free])
            for damping in 2.**-np.arange(9):
                trial = state.copy()
                trial[free] += damping * update
                trial_residual, _ = self.residual(trial, operator, load, tangent=False)
                if np.linalg.norm(trial_residual[free], ord=np.inf) < error:
                    state = trial
                    break
            else:
                raise RuntimeError(f"Newton step stalled with residual {error:.3e}")
        raise RuntimeError(f"Newton reached its iteration limit with residual {history[-1]:.3e}")

    def velocity_boundary(self, nodes, velocity, pressure_reference=None):
        values = velocity(self.velocity_points[nodes])
        fixed = np.r_[nodes, self.n_u + nodes]
        prescribed = np.r_[values[:, 0], values[:, 1]]
        if pressure_reference is not None:
            node, pressure = pressure_reference
            fixed = np.r_[fixed, 2 * self.n_u + node]
            prescribed = np.r_[prescribed, pressure]
        return fixed, prescribed

    def flux_functional(self, edges, normal):
        """Return coefficients for the signed normal flux through tagged edges."""
        lengths = np.linalg.norm(self.points[edges[:, 1]] - self.points[edges[:, 0]], axis=1)
        vector = np.zeros(self.size)
        for i, component in enumerate(normal):
            np.add.at(vector, edges + i * self.n_u, component * lengths[:, None] * np.array([1., 1., 4.]) / 6)
        return vector

    def drag_adjoint(self, state, response_gradient, resistance_gradient):
        """Implicit derivatives for designs entering through quadrature-point drag."""
        free = state.free
        adjoints = np.zeros((self.size, len(response_gradient)))
        factor = splu(state.tangent[free][:, free].tocsc())
        adjoints[free] = factor.solve(np.asarray(response_gradient)[:, free].T, trans="T")
        u, _ = self.velocity_fields(state.values)
        product = np.stack([np.sum(u * self.velocity_fields(adj)[0], axis=-1) for adj in adjoints.T])
        return -np.einsum("req,eq,eqd->rd", product, self.weights, resistance_gradient)
