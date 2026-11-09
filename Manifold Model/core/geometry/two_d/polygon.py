"""An explicitly supplied, straight-sided manifold outline.

No header shape, branch spacing, or dimensions are invented here. The caller
supplies the boundary vertices and labels every edge. Concave outlines are
supported; holes and curved boundaries are not implemented in this first step.
"""

from dataclasses import dataclass
import math

from ..boundary_tags import BoundaryKind, BoundaryTag

Point2D = tuple[float, float]


def _cross(a, b, c):
    """Signed twice-area of triangle a,b,c."""
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _on_segment(a, b, point):
    return (
        _cross(a, b, point) == 0
        and min(a[0], b[0]) <= point[0] <= max(a[0], b[0])
        and min(a[1], b[1]) <= point[1] <= max(a[1], b[1])
    )


def _intersect(a, b, c, d):
    """Detect intersections, including touching and overlapping segments."""
    ab_c, ab_d = _cross(a, b, c), _cross(a, b, d)
    cd_a, cd_b = _cross(c, d, a), _cross(c, d, b)
    if ((ab_c < 0 < ab_d or ab_d < 0 < ab_c)
            and (cd_a < 0 < cd_b or cd_b < 0 < cd_a)):
        return True
    return any((
        _on_segment(a, b, c), _on_segment(a, b, d),
        _on_segment(c, d, a), _on_segment(c, d, b),
    ))


@dataclass(frozen=True)
class PolygonGeometry2D:
    """Counterclockwise boundary vertices in meters, without a repeated end.

    edge_tags[i] labels the edge from vertices[i] to vertices[(i+1) % n].
    Several edges may share a tag. Each inlet/outlet must be one connected
    boundary arc. Wall tags may group multiple disconnected wall segments.
    """

    vertices: tuple[Point2D, ...]
    edge_tags: tuple[BoundaryTag, ...]

    def __post_init__(self):
        # Convert supplied lists to immutable tuples to prevent later edits
        # from silently changing an already configured analysis.
        vertices = tuple(tuple(point) for point in self.vertices)
        tags = tuple(self.edge_tags)
        object.__setattr__(self, "vertices", vertices)
        object.__setattr__(self, "edge_tags", tags)
        count = len(vertices)
        if count < 3 or len(tags) != count:
            raise ValueError("Supply at least three vertices and one tag per edge.")
        if any(len(p) != 2 or not all(math.isfinite(v) for v in p)
               for p in vertices):
            raise ValueError("Vertices must contain finite (x, y) coordinates.")
        if len(set(vertices)) != count:
            raise ValueError("Boundary vertices must be distinct; do not repeat the first.")
        if not all(isinstance(tag, BoundaryTag) for tag in tags):
            raise ValueError("Every edge requires a BoundaryTag.")

        kinds = {}
        for tag in tags:
            if tag.name in kinds and kinds[tag.name] != tag.kind:
                raise ValueError("A boundary name cannot have multiple kinds.")
            kinds[tag.name] = tag.kind

        for i in range(count):
            a, b = vertices[i], vertices[(i + 1) % count]
            # Reject adjacent edges that double back and overlap.
            c = vertices[(i + 2) % count]
            if _on_segment(a, b, c) or _on_segment(b, c, a):
                raise ValueError("Adjacent boundary edges must not overlap.")
            for j in range(i + 1, count):
                if j == i + 1 or (i == 0 and j == count - 1):
                    continue
                if _intersect(a, b, vertices[j], vertices[(j + 1) % count]):
                    raise ValueError("The boundary must not cross or touch itself.")

        if not math.isfinite(self.area_m2) or self.area_m2 <= 0:
            raise ValueError("Boundary vertices must enclose positive counterclockwise area.")

        for name, kind in kinds.items():
            if kind == BoundaryKind.WALL:
                continue
            starts = sum(
                tags[i].name == name and tags[i - 1].name != name
                for i in range(count)
            )
            if starts != 1:
                raise ValueError("Each port must form one connected boundary arc.")

    @property
    def dimension(self):
        return 2

    @property
    def area_m2(self):
        origin = self.vertices[0]
        return 0.5 * math.fsum(
            _cross(origin, self.vertices[i], self.vertices[i + 1])
            for i in range(1, len(self.vertices) - 1)
        )

    def boundary_names(self, kind):
        """Return distinct names in boundary order, with no port-count limit."""
        return tuple(dict.fromkeys(tag.name for tag in self.edge_tags if tag.kind == kind))

    def edge_geometry(self, edge_index):
        """Return start, end, length, and outward unit normal for one edge."""
        if not 0 <= edge_index < len(self.vertices):
            raise IndexError("Boundary edge index is out of range.")
        start = self.vertices[edge_index]
        end = self.vertices[(edge_index + 1) % len(self.vertices)]
        dx, dy = end[0] - start[0], end[1] - start[1]
        length = math.hypot(dx, dy)
        # For a counterclockwise boundary the fluid lies to the left.
        normal = (dy / length, -dx / length)
        return start, end, length, normal
