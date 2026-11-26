"""Common laminar state and discrete-adjoint evaluations for the two models."""

from dataclasses import dataclass
from time import perf_counter

import numpy as np
from scipy.sparse.linalg import splu

from ..project43_power import port_power


HEIGHT_M = 0.004
VISCOSITY = 0.001002
FLUID_DENSITY = 998.2
REYNOLDS = 50.0
VELOCITY_SCALE = REYNOLDS * VISCOSITY / (FLUID_DENSITY * HEIGHT_M)
POWER_SCALE = VISCOSITY * VELOCITY_SCALE**2
RESISTANCE_SCALE = HEIGHT_M**2 / VISCOSITY


@dataclass(frozen=True)
class Stage:
    q: float
    beta: float
    alpha_max: float


@dataclass(frozen=True)
class Evaluation:
    objective: float
    gradient: np.ndarray
    constraints: np.ndarray
    jacobian: np.ndarray
    power_w_per_m: float
    side_fraction: float
    area_mm2: float
    flux: np.ndarray
    mass_imbalance: float
    nonlinear_residual: float
    gray_measure: float


class DesignTransfer:
    """Piecewise-constant cell values and the transpose of their FEM transfer."""

    def __init__(self, cell_index, size):
        self.design_elements = np.asarray(cell_index) >= 0
        self.cells = np.asarray(cell_index)[self.design_elements]
        self.size = size

    def forward(self, dimensional_resistance):
        element_values = np.zeros(len(self.design_elements))
        element_values[self.design_elements] = RESISTANCE_SCALE * dimensional_resistance[self.cells]
        return element_values[:, None]

    def transpose(self, quadrature_cotangent):
        element_cotangent = np.sum(quadrature_cotangent, axis=1)
        return RESISTANCE_SCALE * np.bincount(
            self.cells, weights=element_cotangent[self.design_elements], minlength=self.size
        )


def response_values(domain, state, *, differentiate=True):
    """Return dimensionless port fluxes, physical power, and state gradients."""
    flux = domain.flux @ state
    inlet, main, side = flux
    if not np.all(np.isfinite(flux)) or inlet <= 0:
        raise RuntimeError(f"Invalid signed port fluxes: {flux}")
    imbalance = abs(main + side - inlet) / inlet
    if abs(inlet - 1) > 1e-9 or imbalance > 1e-9:
        raise RuntimeError(f"Port flux balance failed: {flux}, relative imbalance {imbalance:.3e}")
    fraction = side / inlet
    power, gradient = port_power(domain.fem, domain.ports, state, REYNOLDS, differentiate=differentiate)
    state_gradients = (np.stack((POWER_SCALE * gradient, (domain.flux[2] - fraction * domain.flux[0]) / inlet))
                       if differentiate else None)
    return flux, imbalance, POWER_SCALE * sum(power), fraction, state_gradients


class Evaluator:
    """Evaluate one model and continuation stage with a fixed power reference."""

    def __init__(self, domain, pipeline, stage, power_reference, *, reduced=None, initial_state=None):
        if not np.isfinite(power_reference) or power_reference <= 0:
            raise ValueError("Baseline port power must be finite and positive")
        self.domain = domain
        self.pipeline = pipeline
        self.stage = stage
        self.power_reference = power_reference
        self.reduced = reduced
        self.transfer = DesignTransfer(domain.cell_index, pipeline.grid.size)
        self.load = domain.fem.load()
        self.state = initial_state
        self.last_design = None
        self.last_evaluation = None
        self.forward_attempts = 0
        self.forward_solves = 0
        self.newton_iterations = 0
        self.adjoint_rhs = 0
        self.forward_seconds = 0.
        self.adjoint_seconds = 0.

    def raw(self, design):
        if self.reduced is not None:
            return self.reduced.raw(design)
        return design

    def density(self, design):
        return self.pipeline.evaluate(self.raw(design), q=self.stage.q, beta=self.stage.beta,
                                      alpha_max=self.stage.alpha_max)

    def solve(self, resistance, initial):
        self.forward_attempts += 1
        start = perf_counter()
        try:
            state = self.domain.fem.solve(self.domain.fixed, self.domain.prescribed, self.load,
                                          resistance=resistance, initial=initial)
        finally:
            self.forward_seconds += perf_counter() - start
        self.forward_solves += 1
        self.newton_iterations += len(state.residual_history) - 1
        return state

    def design_constraints(self, density, fraction, fraction_gradient):
        area_gradient = self.pipeline.area_gradient(density)
        if self.reduced is not None:
            area_gradient = self.reduced.pullback(area_gradient)
        return (
            np.array([(0.499 - fraction) / 0.001, (fraction - 0.501) / 0.001, density.area_mm2 - 83]),
            np.stack((-fraction_gradient / 0.001, fraction_gradient / 0.001, area_gradient)),
        )

    def evaluate(self, design):
        design = np.asarray(design, dtype=float)
        if self.last_design is not None and np.array_equal(design, self.last_design):
            return self.last_evaluation

        density = self.density(design)
        resistance = self.transfer.forward(density.resistance)
        fem = self.domain.fem
        state = self.solve(resistance, None if self.state is None else self.state.values)
        flux, imbalance, power, fraction, state_gradients = response_values(self.domain, state.values)
        adjoint_start = perf_counter()
        free = state.free
        tangent = splu(state.tangent[free][:, free].tocsc())
        adjoints = np.zeros((fem.size, 2))
        adjoints[free] = tangent.solve(state_gradients[:, free].T, trans="T")
        self.adjoint_rhs += 2
        velocity = fem.velocity_fields(state.values)[0]
        raw_gradients = []
        for adjoint in adjoints.T:
            adjoint_velocity = fem.velocity_fields(adjoint)[0]
            cotangent = -fem.weights * np.sum(velocity * adjoint_velocity, axis=2)
            raw_gradients.append(self.pipeline.pullback(density, self.transfer.transpose(cotangent)))
        if self.reduced is not None:
            raw_gradients = [self.reduced.pullback(gradient) for gradient in raw_gradients]
        self.adjoint_seconds += perf_counter() - adjoint_start
        power_gradient, fraction_gradient = raw_gradients
        constraints, jacobian = self.design_constraints(density, fraction, fraction_gradient)
        gray = float(4 * np.mean(density.projected * (1 - density.projected)))
        result = Evaluation(
            power / self.power_reference,
            power_gradient / self.power_reference,
            constraints,
            jacobian,
            power, fraction, density.area_mm2, flux, imbalance, state.residual_history[-1], gray,
        )
        self.state = state
        self.last_design = design.copy()
        self.last_evaluation = result
        return result


def baseline_reference(domain, pipeline, alpha_max):
    """Return common P0 and the converged baseline state using inertia continuation."""
    stage = Stage(0.01, 0.0, alpha_max)
    evaluator = Evaluator(domain, pipeline, stage, 1.0)
    density = evaluator.density(pipeline.grid.base_raw)
    resistance = evaluator.transfer.forward(density.resistance)
    fem = domain.fem
    state = None
    try:
        for inertia in (0.0, REYNOLDS / 2, REYNOLDS):
            fem.inertia = inertia
            state = evaluator.solve(resistance, None if state is None else state.values)
    finally:
        fem.inertia = REYNOLDS
    power = response_values(domain, state.values)[2]
    if not np.isfinite(power) or power <= 0:
        raise RuntimeError(f"Nonpositive baseline port power: {power} W/m")
    return power, state
