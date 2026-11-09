"""Explicit inlet and outlet data; no default pressures or flow rates.

Pressure outlet values specify the pressure load in sigma*n = -p_load*n.
The linear-resistance option follows the relationship in Outline/Analysis.md:
p_load = p_reference + resistance * Q_out, with outward-positive Q_out.
"""

from collections.abc import Callable
from dataclasses import dataclass
import math


@dataclass(frozen=True)
class InletVelocity2D:
    """A supplied depth-averaged vector profile U(x,y), in m/s."""

    boundary_name: str
    profile: Callable[[float, float], tuple[float, float]]

    def __post_init__(self):
        if not isinstance(self.boundary_name, str) or not self.boundary_name.strip():
            raise ValueError("The inlet must name a geometry boundary.")
        if not callable(self.profile):
            raise ValueError("Supply a callable inlet velocity profile.")

    def velocity_at(self, x, y):
        values = tuple(self.profile(x, y))
        if len(values) != 2 or not all(math.isfinite(v) for v in values):
            raise ValueError("The inlet profile must return two finite velocity components.")
        return values


@dataclass(frozen=True)
class PressureOutlet:
    boundary_name: str
    pressure_load_pa: float

    def __post_init__(self):
        if not isinstance(self.boundary_name, str) or not self.boundary_name.strip():
            raise ValueError("The outlet must name a geometry boundary.")
        if not math.isfinite(self.pressure_load_pa):
            raise ValueError("Outlet pressure must be finite.")


@dataclass(frozen=True)
class LinearResistanceOutlet:
    boundary_name: str
    reference_pressure_pa: float
    resistance_pa_s_per_m3: float

    def __post_init__(self):
        if not isinstance(self.boundary_name, str) or not self.boundary_name.strip():
            raise ValueError("The outlet must name a geometry boundary.")
        if not math.isfinite(self.reference_pressure_pa):
            raise ValueError("Reference pressure must be finite.")
        if (not math.isfinite(self.resistance_pa_s_per_m3)
                or self.resistance_pa_s_per_m3 < 0):
            raise ValueError("Downstream resistance must be finite and nonnegative.")


OutletCondition = PressureOutlet | LinearResistanceOutlet
