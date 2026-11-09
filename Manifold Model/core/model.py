"""Connect explicit geometry, fluid, and boundary inputs for one fixed model."""

from dataclasses import dataclass
import math

from .boundaries.conditions import (
    InletVelocity2D, LinearResistanceOutlet, OutletCondition, PressureOutlet,
)
from .geometry.boundary_tags import BoundaryKind
from .geometry.two_d.polygon import PolygonGeometry2D
from .physics.material import FluidProperties


@dataclass(frozen=True)
class DepthAveragedModel2D:
    """Steady incompressible laminar analysis of a constant-depth manifold.

    This selects the depth-averaged model described in Outline/Analysis.md.
    Walls are stationary with no slip. Geometry stays fixed. Exactly one
    inlet and any positive number of outlets are supported at this stage.
    This object describes inputs; it does not mesh or solve the problem.
    """

    geometry: PolygonGeometry2D
    depth_m: float
    fluid: FluidProperties
    inlet: InletVelocity2D
    outlets: tuple[OutletCondition, ...]

    def __post_init__(self):
        object.__setattr__(self, "outlets", tuple(self.outlets))
        if not isinstance(self.geometry, PolygonGeometry2D):
            raise ValueError("Supply a PolygonGeometry2D.")
        if not isinstance(self.fluid, FluidProperties):
            raise ValueError("Supply constant FluidProperties.")
        if not isinstance(self.inlet, InletVelocity2D):
            raise ValueError("Supply an InletVelocity2D condition.")
        if not math.isfinite(self.depth_m) or self.depth_m <= 0:
            raise ValueError("Manifold depth must be finite and positive.")
        inlets = self.geometry.boundary_names(BoundaryKind.INLET)
        outlets = self.geometry.boundary_names(BoundaryKind.OUTLET)
        if inlets != (self.inlet.boundary_name,):
            raise ValueError("Label exactly one inlet and supply its velocity condition.")
        if not outlets or not self.geometry.boundary_names(BoundaryKind.WALL):
            raise ValueError("The manifold needs outlet and wall boundaries.")
        if not all(isinstance(o, (PressureOutlet, LinearResistanceOutlet))
                   for o in self.outlets):
            raise ValueError("Supply a supported condition for each outlet.")
        names = tuple(o.boundary_name for o in self.outlets)
        if len(set(names)) != len(names) or set(names) != set(outlets):
            raise ValueError("Supply exactly one condition for every labeled outlet.")

    @property
    def outlet_count(self):
        return len(self.outlets)

    def target_outlet_flow(self, inlet_flow_m3_per_s):
        """Equal-flow target for reporting only; this does not enforce flow."""
        if not math.isfinite(inlet_flow_m3_per_s) or inlet_flow_m3_per_s <= 0:
            raise ValueError("Supply positive incoming volumetric flow.")
        return inlet_flow_m3_per_s / self.outlet_count
