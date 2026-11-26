"""Flow and power evaluation for the prescribed wall-displacement experiment."""

from copy import copy
from dataclasses import replace

import numpy as np

from .evaluator import Evaluator


PDE_MODES = ("navier-stokes", "stokes")


class WallDisplacementEvaluatorV1(Evaluator):
    """Reuse the original FEM and two adjoints with a flow-only constraint set."""

    def __init__(self, domain, pipeline, stage, power_reference, *, initial_state=None, pde="navier-stokes"):
        if pde not in PDE_MODES:
            raise ValueError(f"Unsupported wall-displacement PDE: {pde}")
        if pde == "stokes":
            fem = copy(domain.fem)
            fem.inertia = 0.
            domain = replace(domain, fem=fem)
        super().__init__(domain, pipeline, stage, power_reference, initial_state=initial_state)
        self.pde = pde

    def design_constraints(self, density, fraction, fraction_gradient):
        sensitivity = fraction_gradient / .001
        return (
            np.array([(.10 - fraction) / .001, (fraction - .90) / .001]),
            np.stack((-sensitivity, sensitivity)),
        )
