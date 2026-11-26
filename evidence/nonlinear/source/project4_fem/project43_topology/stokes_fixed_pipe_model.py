"""Fixed-pipe Stokes variant with the same physical port responses."""

from .fixed_pipe_pilot_model import SPACING_MM, FixedPipeEvaluator, build_fixed_pipe


def build_stokes_fixed_pipe(max_spacing_mm=SPACING_MM, density_spacing_mm=SPACING_MM):
    """Keep the pilot geometry and boundary data and remove state advection."""
    domain, pipeline = build_fixed_pipe(max_spacing_mm, density_spacing_mm)
    domain.fem.inertia = 0.
    return domain, pipeline


class StokesFixedPipeEvaluator(FixedPipeEvaluator):
    """Use the existing affine Taylor–Hood residual and discrete adjoint."""

    def __init__(self, domain, pipeline, stage, power_reference, side_reference, *, initial_state=None):
        if domain.fem.inertia != 0.:
            raise ValueError("A Stokes evaluator requires a zero-inertia FEM domain")
        super().__init__(domain, pipeline, stage, power_reference, side_reference, initial_state=initial_state)
