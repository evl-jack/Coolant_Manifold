"""P2/P1 assembly for constant-depth, incompressible 2D verification cases."""

from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import spsolve


def line_rule(order):
    """Map Gauss–Legendre points and weights to the unit interval."""
    points, weights = leggauss(order)
    return (1 + points) / 2, weights / 2


def triangle_rule(order):
    """Map tensor-product quadrature onto the reference triangle."""
    t, w = line_rule(order)
    a, b = np.meshgrid(t, t, indexing="ij")
    return np.column_stack((a.ravel(), ((1 - a) * b).ravel())), ((1 - a) * np.outer(w, w)).ravel()


def basis(points):
    """Evaluate P1 and P2 bases and reference gradients at (xi, eta)."""
    xi, eta = points.T
    l = np.column_stack((1 - xi - eta, xi, eta))
    dl = np.array([[-1, -1], [1, 0], [0, 1]])
    first, second = [0, 1, 2], [1, 2, 0]
    p2 = np.column_stack((l * (2*l - 1), 4*l[:, first]*l[:, second]))
    dp2 = np.concatenate(((4*l - 1)[..., None]*dl,
                          4*(l[:, first, None]*dl[second] + l[:, second, None]*dl[first])), axis=1)
    return l, p2, dp2


@dataclass
class RectangleGrid:
    """Build shared nodes and element connectivity for both trial spaces."""
    nx: int
    ny: int

    def __post_init__(self):
        # Store horizontal velocity, vertical velocity, pressure, then the gauge.
        self.p = self.unit_nodes(self.nx, self.ny)
        self.v = self.unit_nodes(2*self.nx, 2*self.ny)
        self.p_tri = self.elements(1)
        self.v_tri = self.elements(2)
        self.nv, self.np = len(self.v), len(self.p)
        self.ids = np.column_stack((self.v_tri, self.nv + self.v_tri, 2*self.nv + self.p_tri))
        self.size = 2*self.nv + self.np + 1
        self.kind = np.tile([0, 1], self.nx*self.ny)

    @staticmethod
    def unit_nodes(nx, ny):
        x, y = np.meshgrid(np.linspace(0, 1, nx+1), np.linspace(0, 1, ny+1))
        return np.column_stack((x.ravel(), y.ravel()))

    def elements(self, degree):
        # Split each rectangle diagonally and add edge midpoints for P2.
        width = degree*self.nx + 1
        x, y = np.meshgrid(degree*np.arange(self.nx), degree*np.arange(self.ny))
        origin = (y*width + x).ravel()
        dx, dy = degree, degree*width
        offsets = [[0, dx, dx+dy], [0, dx+dy, dy]]
        if degree == 2:
            offsets[0] += [1, 2+width, 1+width]
            offsets[1] += [1+width, 1+2*width, width]
        return (origin[:, None, None] + np.asarray(offsets)).reshape(-1, len(offsets[0]))


@dataclass
class RectangleMesh:
    """Map the normalized grid into physical coordinates, with optional rotation."""
    grid: RectangleGrid
    length: float
    height: float
    rotated: bool = False

    def __post_init__(self):
        # Each rectangle contributes two affine triangle mappings.
        self.rotation = np.array([[0, -1], [1, 0]]) if self.rotated else np.eye(2)
        dimensions = np.array([self.length, self.height])
        self.p = (self.grid.p*dimensions) @ self.rotation.T
        self.v = (self.grid.v*dimensions) @ self.rotation.T
        dx, dy = dimensions / [self.grid.nx, self.grid.ny]
        self.jacobian = self.rotation @ np.array([[[dx, dx], [0, dy]], [[dx, 0], [dy, dy]]])
        self.inverse = np.linalg.inv(self.jacobian)
        self.det = dx*dy
        self.area = self.length*self.height
        self.h_max = np.hypot(dx, dy)

    def physical_points(self, reference_points):
        """Return the physical quadrature points for every element."""
        origins = self.p[self.grid.p_tri[:, 0]]
        offsets = np.einsum("kij,qj->kqi", self.jacobian, reference_points)
        return origins[:, None, :] + offsets[self.grid.kind]


def element_matrix(grid, local):
    """Scatter local element matrices into a global sparse matrix."""
    rows = np.broadcast_to(grid.ids[:, :, None], local.shape)
    cols = np.broadcast_to(grid.ids[:, None, :], local.shape)
    return coo_matrix((local.ravel(), (rows.ravel(), cols.ravel())),
                      shape=(grid.size, grid.size)).tocsr()


def scatter(grid, local):
    """Sum element vector entries at their shared global degrees of freedom."""
    return np.bincount(grid.ids.ravel(), local.ravel(), minlength=grid.size)


def strain_and_divergence(grad):
    """Build sqrt(2) normal strains, engineering shear, and divergence."""
    dx, dy = grad[..., 0], grad[..., 1]
    zero = np.zeros_like(dx)
    join = lambda a, b: np.concatenate((a, b), axis=-1)
    strain = np.stack((join(np.sqrt(2)*dx, zero), join(zero, np.sqrt(2)*dy), join(dy, dx)), axis=-2)
    return strain, join(dx, dy)


class Assembly:
    """Build weak-form operators on a fixed mesh using cached basis evaluations."""
    def __init__(self, mesh, order=6):
        self.mesh = mesh
        self.grid = mesh.grid
        self.points, self.weights = triangle_rule(order)
        self.p, self.n, grad_ref = basis(self.points)
        self.grad = grad_ref[None, ...] @ mesh.inverse[:, None, ...]
        self.strain, self.div = strain_and_divergence(self.grad)
        self.mass = np.einsum("q,qi,qj->ij", self.weights, self.n, self.n)

    def linear(self, mu, drag, force, direction=None):
        """Assemble a system, or its derivative for an axis-aligned uniform strip.

        The optional direction contains derivatives of mu and drag, logarithmic
        changes of L and H, and the partial of the volume force at fixed X,Y.
        """
        m, g, w = self.mesh, self.grid, self.weights
        local = np.zeros((2, 15, 15))
        # Symmetric-gradient viscosity and the pressure/divergence coupling.
        visc = np.einsum("q,kqai,kqaj->kij", w, self.strain, self.strain)
        pressure = -np.einsum("q,qi,kqj->kij", w, self.p, self.div)
        determinant_rate = 0.0
        if direction is None:
            viscosity_block = mu*visc
            drag_factor = drag
            values = force(*np.moveaxis(m.physical_points(self.points), -1, 0))
        else:
            # Differentiate both the integration measure and physical gradients
            # when the rectangular geometry changes at fixed normalized nodes.
            dmu, ddrag, log_l, log_h, dforce = direction
            determinant_rate = log_l + log_h
            dstrain, ddiv = strain_and_divergence(-self.grad * [log_l, log_h])
            cross = np.einsum("q,kqai,kqaj->kij", w, dstrain, self.strain)
            viscosity_block = (dmu + mu*determinant_rate)*visc + mu*(cross + cross.transpose(0, 2, 1))
            drag_factor = ddrag + drag*determinant_rate
            pressure = pressure*determinant_rate - np.einsum("q,qi,kqj->kij", w, self.p, ddiv)
            x, y = np.moveaxis(m.physical_points(self.points), -1, 0)
            values = np.asarray(dforce(x, y)) + determinant_rate*np.asarray(force(x, y))

        local[:, :12, :12] = viscosity_block
        local[:, :6, :6] += drag_factor*self.mass
        local[:, 6:12, 6:12] += drag_factor*self.mass
        local[:, 12:, :12] = pressure
        local[:, :12, 12:] = pressure.transpose(0, 2, 1)
        matrix = element_matrix(g, m.det*local[g.kind])
        if direction is None:
            # Append the normalized mean-pressure constraint and its multiplier.
            gauge = np.bincount(g.p_tri.ravel(), np.tile(np.full(3, m.det/(6*m.area)), len(g.ids)), minlength=g.np)
            rows = np.arange(2*g.nv, 2*g.nv + g.np)
            last = np.full(g.np, g.size - 1)
            matrix += coo_matrix((np.tile(gauge, 2), (np.r_[rows, last], np.r_[last, rows])), shape=matrix.shape).tocsr()
        rhs = np.zeros((len(g.ids), 15))
        rhs[:, :12] = (m.det*np.einsum("q,qi,aeq->eai", w, self.n, values)).reshape(-1, 12)
        return matrix, scatter(g, rhs)

    def velocity_on_rule(self, full):
        """Interpolate velocity and its physical gradient at quadrature points."""
        g = self.grid
        uv = full[:2*g.nv].reshape(2, -1).T[g.v_tri]
        value = np.einsum("qi,eia->eqa", self.n, uv)
        gradient = np.einsum("eqib,eia->eqab", self.grad[g.kind], uv)
        return value, gradient

    def convection(self, full, density, tangent):
        """Assemble the inertial residual and optionally its Newton derivative."""
        g, m = self.grid, self.mesh
        velocity, gradient = self.velocity_on_rule(full)
        w = density*m.det*self.weights
        advected = np.einsum("eqab,eqb->eqa", gradient, velocity)
        local = np.zeros((len(g.ids), 15))
        local[:, :12] = np.einsum("q,qi,eqa->eai", w, self.n, advected).reshape(-1, 12)
        if not tangent:
            return scatter(g, local)
        matrix = np.zeros((len(g.ids), 15, 15))
        # The derivative varies both the advecting velocity and the advected field.
        advector = np.einsum("q,qi,eqab,qj->eaibj", w, self.n, gradient, self.n)
        advected_derivative = np.einsum("q,qi,eqb,eqjb->eij", w, self.n, velocity, self.grad[g.kind])
        advector[:, 0, :, 0, :] += advected_derivative
        advector[:, 1, :, 1, :] += advected_derivative
        matrix[:, :12, :12] = advector.reshape(-1, 12, 12)
        return scatter(g, local), element_matrix(g, matrix)


def vector_constant(x, y, first=0.0, second=0.0):
    """Broadcast a constant body force over the element quadrature points."""
    return np.broadcast_to(np.array([first, second])[:, None, None], (2, *x.shape))


class Problem:
    """Reduce the equations to free unknowns after applying boundary data."""
    def __init__(self, assembly, mu, drag, force, *, periodic=False, boundary=None, density=0):
        self.assembly, self.density = assembly, density
        g = assembly.grid
        a, b = assembly.linear(mu, drag, force)
        reps = np.arange(g.size)
        x, y = g.v.T
        fixed = (y == 0) | (y == 1)
        if periodic:
            # Identify matching left and right degrees of freedom, including pressure.
            for offset, points in [(0, g.v), (g.nv, g.v), (2*g.nv, g.p)]:
                left, right = points[:, 0] == 0, points[:, 0] == 1
                reps[offset + np.flatnonzero(right)] = offset + np.flatnonzero(left)
        else:
            fixed |= (x == 0) | (x == 1)
        free = np.ones(g.size, dtype=bool)
        free[:2*g.nv] = ~np.tile(fixed, 2)
        lift = np.zeros(g.size)
        if boundary is not None:
            lift[np.flatnonzero(~free)] = np.asarray(boundary(x[fixed], y[fixed])).ravel()
        _, cols = np.unique(reps[free], return_inverse=True)
        # Recover the full field as prolong @ state + lift. Its transpose
        # assembles the residual on the corresponding free equations.
        self.prolong = coo_matrix((np.ones(len(cols)), (np.flatnonzero(free), cols)),
                                  shape=(g.size, int(cols.max())+1)).tocsr()
        self.transpose = self.prolong.T.tocsr()
        self.lift = lift
        self.matrix = (self.transpose @ a @ self.prolong).tocsr()
        self.rhs = self.transpose @ (b - a @ lift)

    def recover(self, state):
        """Restore prescribed boundary values and periodic copies."""
        return self.prolong @ state + self.lift

    def residual(self, state):
        """Evaluate the reduced momentum, continuity, and gauge residual."""
        nonlinear = self.assembly.convection(self.recover(state), self.density, False)
        return self.matrix @ state - self.rhs + self.transpose @ nonlinear

    def jacobian(self, state):
        """Differentiate the reduced residual with respect to the free state."""
        _, nonlinear = self.assembly.convection(self.recover(state), self.density, True)
        return self.matrix + self.transpose @ nonlinear @ self.prolong


def backward_error(matrix, state, rhs):
    """Measure the residual relative to the scale of the linear system."""
    denom = (np.asarray(abs(matrix).sum(axis=1)).max()*np.linalg.norm(state, np.inf)
             + np.linalg.norm(rhs, np.inf))
    residual = np.linalg.norm(matrix @ state - rhs, np.inf)
    return float(residual / denom if denom else residual)


def direct(matrix, rhs):
    """Equilibrate and solve a sparse linear system, then report backward error."""
    scale = diags(1/np.sqrt(abs(matrix).max(axis=1).toarray().ravel()))
    state = scale @ spsolve((scale @ matrix @ scale).tocsc(), scale @ rhs)
    if not np.isfinite(state).all():
        raise ArithmeticError("Sparse solve returned a nonfinite state")
    return state, backward_error(matrix, state, rhs)


def newton(problem, atol=1e-15, rtol=1e-12, max_steps=25, max_backtracks=20,
           *, initial_state=None):
    """Solve the reduced nonlinear equations with backtracking Newton steps.

    The default starts from the Stokes solution. Verification can supply an
    initial reduced state and request just one accepted step.
    """
    if initial_state is None:
        state, initial_error = direct(problem.matrix, problem.rhs)
    else:
        state = np.array(initial_state, dtype=float, copy=True)
        if state.shape != problem.rhs.shape or not np.isfinite(state).all():
            raise ValueError("Initial Newton state must be a finite reduced state.")
        initial_error = None
    residual = problem.residual(state)
    history = [float(np.linalg.norm(residual))]
    threshold = atol + rtol*max(history[0], np.linalg.norm(problem.rhs))
    steps, errors = [], []
    while history[-1] > threshold and len(steps) < max_steps:
        correction, error = direct(problem.jacobian(state), -residual)
        # Try the full Newton update first and halve it if the residual
        # fails the sufficient-decrease test.
        for backtrack in range(max_backtracks+1):
            step = 0.5**backtrack
            trial = state + step*correction
            trial_residual = problem.residual(trial)
            norm = np.linalg.norm(trial_residual)
            if norm <= (1 - 1e-4*step)*history[-1]:
                state, residual = trial, trial_residual
                steps.append(step)
                errors.append(error)
                history.append(float(norm))
                break
        else:
            raise ArithmeticError("Newton backtracking exhausted")
    return state, {"converged": bool(history[-1] <= threshold), "history": history, "steps": steps,
                            "atol": atol, "rtol": rtol, "threshold": threshold,
                            "initial_backward_error": initial_error, "step_backward_errors": errors}
