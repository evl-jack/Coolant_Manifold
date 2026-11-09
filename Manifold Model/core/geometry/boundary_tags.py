"""Boundary names connect geometry to conditions without hard-coded port counts."""

from dataclasses import dataclass
from enum import Enum


class BoundaryKind(str, Enum):
    INLET = "inlet"
    OUTLET = "outlet"
    WALL = "wall"


@dataclass(frozen=True)
class BoundaryTag:
    """One physical boundary name, such as a port or a group of walls."""

    name: str
    kind: BoundaryKind

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("A boundary needs a nonempty name.")
        if not isinstance(self.kind, BoundaryKind):
            raise ValueError("Boundary kind must be a BoundaryKind value.")
