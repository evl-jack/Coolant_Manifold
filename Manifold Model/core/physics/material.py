"""Constant physical properties for the specified incompressible fluid."""

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class FluidProperties:
    density_kg_per_m3: float
    dynamic_viscosity_pa_s: float

    def __post_init__(self):
        for value in (self.density_kg_per_m3, self.dynamic_viscosity_pa_s):
            if not math.isfinite(value) or value <= 0:
                raise ValueError("Fluid density and viscosity must be finite and positive.")
