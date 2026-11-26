"""Differentiable paired fillet design and its reduced FEM responses."""

import numpy as np
from scipy.interpolate import BSpline
from scipy.sparse import coo_matrix

from .fem import TaylorHood
from .junction import boundary_data, junction_mesh


CONTROLS = 20
BASELINE = np.full(CONTROLS, .5)
ALPHA = 20_000.
TOLERANCE = .005
RESPONSES = ("dissipation", "side_share", "pressure_loss", "solid_area_mm2")


class FilletMap:
    def __init__(self, geometry):
        self.main, self.outlet, self.stem, self.radius = geometry
        self.root = 4. + self.main
        knots = np.r_[np.zeros(4), np.linspace(0., 1., CONTROLS-2)[1:-1], np.ones(4)]
        self.spline = BSpline(knots, np.eye(CONTROLS), 3)

    def parameters(self, y):
        h = np.clip((y-self.root)/self.radius, 0., 1.)
        s = np.arccos(1-h)/(np.pi/2)
        clearance = self.radius*(1-np.sin(.5*np.pi*s))
        raw = .45*np.sin(np.pi*s)**2
        a = np.divide(raw*clearance, raw+clearance, out=np.zeros_like(raw), where=raw+clearance > 0)
        return s, clearance, a

    def offset(self, y, design=None):
        s, _, a = self.parameters(y)
        if design is None:
            return 1.14*a
        return a*(self.spline(s)@(2*design-1))

    def prepare(self, points):
        x = points[..., 0]*(2*self.main)-1
        y = points[..., 1]*(2*self.main)+4-self.main
        s, clearance, a = self.parameters(y)
        self.fillet = np.argwhere((self.root < y) & (y < self.root+self.radius) & (a > 0))
        self.index = tuple(self.fillet.T)
        self.shape = y.shape
        distance = abs(x-10)-self.stem-clearance
        self.edit = (self.root < y) & (y < self.root+self.radius) & (abs(distance) <= 1.14*a)
        transition = .28*a[self.index]
        self.origin = .5-distance[self.index]/transition
        self.design_basis = self.spline(s[self.index])/.14

    def density(self, design):
        z = np.clip(self.origin+self.design_basis@(design-.5), 0., 1.)
        gamma = np.ones(self.shape)
        gamma[self.index] = z*z*(3-2*z)
        slope = (6*z*(1-z))[:, None]*self.design_basis
        return gamma, slope


def pressure_weights(fem, ports):
    """Inlet pressure-mean vector and outlet pressure-times-flux bilinear form."""
    inlet, *outlets = ports
    length = np.linalg.norm(fem.points[inlet[:, 1]]-fem.points[inlet[:, 0]], axis=1)
    reference = np.zeros(fem.size)
    np.add.at(reference, 2*fem.n_u+inlet[:, :2], length[:, None]*np.ones((1, 2))/(2*length.sum()))
    # Exact line integral of the P2 velocity basis times the P1 pressure basis.
    product = np.array([[1/6, 0.], [0., 1/6], [1/3, 1/3]])
    matrix = coo_matrix((fem.size, fem.size)).tocsr()
    for i, port in enumerate(outlets):
        length = np.linalg.norm(fem.points[port[:, 1]]-fem.points[port[:, 0]], axis=1)
        data = length[:, None, None]*product
        u = np.broadcast_to((i*fem.n_u+port)[:, :, None], data.shape)
        p = np.broadcast_to((2*fem.n_u+port[:, :2])[:, None, :], data.shape)
        matrix += coo_matrix((data.ravel(), (u.ravel(), p.ravel())), shape=matrix.shape).tocsr()
    return reference, matrix


class Problem:
    def __init__(self, geometry, reynolds=50., refine=1, alpha=ALPHA, *, mapping=None, transition_rows=6):
        self.geometry, self.reynolds, self.alpha = geometry, reynolds, alpha
        self.map = FilletMap(geometry) if mapping is None else mapping
        self.fem = TaylorHood(*junction_mesh(geometry, refine, self.map.offset, transition_rows), inertia=reynolds)
        self.map.prepare(self.fem.quad_points)
        self.fixed, self.values, self.ports, self.flux = boundary_data(self.fem)
        self.p_in, self.p_out = pressure_weights(self.fem, self.ports)
        self.load = self.fem.load()
        self.viscous = self.fem.viscous_pressure[:2*self.fem.n_u, :2*self.fem.n_u]
        self.area_weights = self.fem.weights*self.map.edit*(2*self.map.main)**2
        self.warm = None
        self.cached_design = None

    def evaluate(self, design, gradients=False):
        fem, idx = self.fem, self.map.index
        if self.cached_design is None or not np.array_equal(design, self.cached_design):
            gamma, slope = self.map.density(design)
            resistance = self.alpha*(1-gamma)**2
            if self.warm is None and self.reynolds:
                fem.inertia = 0.
                try:
                    self.warm = fem.solve(self.fixed, self.values, self.load, resistance=resistance).values
                finally:
                    fem.inertia = self.reynolds
            state = fem.solve(self.fixed, self.values, self.load, resistance=resistance, initial=self.warm)
            self.warm = state.values
            self.cached_design = np.array(design, copy=True)
            self.cache = gamma, slope, resistance, state
        gamma, slope, resistance, state = self.cache
        w = state.values
        u, _ = fem.velocity_fields(w)
        speed2 = np.sum(u*u, axis=-1)
        u_vector = w[:2*fem.n_u]
        viscous = self.viscous @ u_vector
        dissipation = u_vector @ viscous + np.sum(fem.weights*resistance*speed2)
        q_in, q_m, q_s = self.flux @ w
        out_flux = q_m+q_s
        p_weighted = w @ (self.p_out@w)
        response = np.array([dissipation, q_s/q_in, self.p_in@w-p_weighted/out_flux,
                             np.sum(self.area_weights*(1-gamma))])
        if not gradients:
            return response

        partials = np.zeros((3, fem.size))
        partials[0, :2*fem.n_u] = 2*viscous
        partials[0] += 2*fem.scatter_velocity(np.einsum("eq,qa,eqi->eai", fem.weights*resistance, fem.basis, u))
        partials[1] = (self.flux[2]-response[1]*self.flux[0])/q_in
        partials[2] = (self.p_in-(self.p_out+self.p_out.T)@w/out_flux
                       + p_weighted*(self.flux[1]+self.flux[2])/out_flux**2)
        resistance_gradient = np.zeros((*gamma.shape, slope.shape[1]))
        resistance_gradient[idx] = -2*self.alpha*(1-gamma[idx])[:, None]*slope
        reduced = fem.drag_adjoint(state, partials, resistance_gradient)
        reduced[0] += (fem.weights[idx]*speed2[idx])@resistance_gradient[idx]
        return response, np.vstack((reduced, -self.area_weights[idx]@slope))

    def diagnostics(self, design):
        self.evaluate(design)
        gamma, _, _, state = self.cache
        u, _ = self.fem.velocity_fields(state.values)
        q = self.flux @ state.values
        solid = gamma < .05
        weights = self.fem.weights[solid]
        solid_rms = float(np.sqrt(weights @ np.sum(u[solid]**2, axis=-1)/weights.sum())) if weights.size else None
        return dict(mass_imbalance=float(abs(q[0]-q[1]-q[2])/q[0]), outlet_fluxes=q[1:].tolist(),
                    solid_rms_speed=solid_rms, nominal_solid_quadrature_points=int(solid.sum()),
                    residual=state.residual_history[-1])
