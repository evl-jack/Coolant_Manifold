"""A supplied planar outline extruded through a supplied physical depth.

This is geometry only, not a mesh or a fluid solver. Face ordering is outward.
The first geometry family preserves the existing 2D port names and count.
"""

from dataclasses import dataclass
import math

from ..boundary_tags import BoundaryKind, BoundaryTag
from ..two_d.polygon import PolygonGeometry2D


@dataclass(frozen=True)
class ExtrudedGeometry3D:
    outline: PolygonGeometry2D
    depth_m: float
    base_z_m: float
    bottom_wall: BoundaryTag
    top_wall: BoundaryTag

    def __post_init__(self):
        if not isinstance(self.outline, PolygonGeometry2D):
            raise ValueError("Supply the actual planar outline.")
        if not math.isfinite(self.depth_m) or self.depth_m <= 0:
            raise ValueError("Extrusion depth must be finite and positive.")
        if not math.isfinite(self.base_z_m):
            raise ValueError("The base elevation must be finite.")
        top = self.base_z_m + self.depth_m
        if not math.isfinite(top) or top <= self.base_z_m:
            raise ValueError("The extrusion must have a distinct finite top elevation.")
        for tag in (self.bottom_wall, self.top_wall):
            if not isinstance(tag, BoundaryTag) or tag.kind != BoundaryKind.WALL:
                raise ValueError("The explicit top and bottom must be wall boundaries.")
            for side in self.outline.edge_tags:
                if tag.name == side.name and side.kind != BoundaryKind.WALL:
                    raise ValueError("A wall name cannot replace a port name.")

    @property
    def dimension(self):
        return 3

    @property
    def vertices(self):
        """Bottom vertices first, then corresponding top vertices, in meters."""
        return tuple(
            (x, y, z)
            for z in (self.base_z_m, self.base_z_m + self.depth_m)
            for x, y in self.outline.vertices
        )

    @property
    def faces(self):
        """Return (vertex_ids, boundary_tag) for the enclosing polygon faces.

        Top/bottom faces may be concave polygons. They are not triangle fans
        and must be triangulated by a mesher that supports the actual outline.
        """
        n = len(self.outline.vertices)
        faces = [
            (tuple(reversed(range(n))), self.bottom_wall),
            (tuple(range(n, 2 * n)), self.top_wall),
        ]
        for i, tag in enumerate(self.outline.edge_tags):
            j = (i + 1) % n
            faces.append(((i, j, j + n, i + n), tag))
        return tuple(faces)

    @property
    def volume_m3(self):
        return self.outline.area_m2 * self.depth_m

    def boundary_names(self, kind):
        return tuple(dict.fromkeys(tag.name for _, tag in self.faces if tag.kind == kind))

    def vertex_depth_derivative(self):
        """dX/db at fixed outline and fixed bottom elevation."""
        n = len(self.outline.vertices)
        return ((0.0, 0.0, 0.0),) * n + ((0.0, 0.0, 1.0),) * n

    def volume_depth_derivative(self):
        """dV/db at fixed outline."""
        return self.outline.area_m2
