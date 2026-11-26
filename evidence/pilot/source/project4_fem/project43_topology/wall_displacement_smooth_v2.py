"""Single constrained-wall experiment using the existing threshold reconstruction."""

import numpy as np
from scipy.special import logsumexp
from shapely.geometry import LineString, MultiLineString

from .reconstruction import _PORTS_MM, reconstruct
from .wall_displacement_density_v1 import WallDisplacementDensityV1
from .wall_displacement_geometry_v1 import _clipped, editable_geometry, sharp_geometry


STEP = 1e-5
SECT = .25
SAMPLES = 1001
TAU = .001
LIMIT = np.cos(np.pi / 4)
SCALE = .1


class ProtectedDensity(WallDisplacementDensityV1):
    core_boxes = ((0., 0., 20., 3.75), (9.5, 2., 10.5, 7.))


class UnresolvedGeometry(ValueError):
    pass


def wall_chains(geometry):
    if geometry.port_component is None or len(geometry.polygons) != 1 or geometry.polygons[0].interiors:
        raise UnresolvedGeometry("Wall chains require a single port-connected fluid polygon without holes")
    ring = np.asarray(geometry.polygons[0].exterior.coords[:-1])
    endpoints = {tuple(end) for port in _PORTS_MM for end in port}
    tips = [i for i, point in enumerate(ring) if tuple(point) in endpoints]
    if len(tips) != 6:
        raise UnresolvedGeometry("Could not split the wall at all six exact port endpoints")
    chains = []
    for start, end in zip(tips, tips[1:] + [tips[0] + len(ring)]):
        points = ring[np.arange(start, end + 1) % len(ring)]
        terminal = tuple(sorted((tuple(points[0]), tuple(points[-1]))))
        if terminal not in {tuple(sorted(port)) for port in _PORTS_MM}:
            chains.append((terminal, points if tuple(points[0]) == terminal[0] else points[::-1]))
    if len(chains) != 3:
        raise UnresolvedGeometry("Expected exactly three port-to-port wall chains")
    return [points for _, points in sorted(chains, key=lambda chain: chain[0])]


def arc_table(points):
    arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))]
    unique = np.r_[True, np.diff(arc) > 1e-12]
    return arc[unique], points[unique]


def cosines(points, arc, positions):
    sampled = [np.column_stack([np.interp(np.clip(positions + offset, 0., arc[-1]), arc, axis)
                                for axis in points.T]) for offset in (-SECT, 0., SECT)]
    previous, middle, following = sampled
    incoming, outgoing = middle - previous, following - middle
    scale = np.linalg.norm(incoming, axis=1) * np.linalg.norm(outgoing, axis=1)
    if np.any(scale <= 1e-20):
        raise UnresolvedGeometry("An eligible wall sample produced a degenerate secant")
    return np.clip(np.einsum("ij,ij->i", incoming, outgoing) / scale, -1., 1.)


def sample_angles(chains, *, verification=False):
    slopes, spacing = [], []
    for points in chains:
        arc, points = arc_table(points)
        lower, upper = SECT, arc[-1] - SECT
        if upper <= lower:
            raise UnresolvedGeometry("A wall chain has no valid secant interval")
        step = (upper - lower) / (SAMPLES - 1)
        spacing.append(float(step))
        if step > .025 + 1e-12:
            raise UnresolvedGeometry(f"A wall chain has {step:.8g} mm sample spacing, above 0.025 mm")
        samples = np.linspace(lower, upper, SAMPLES)
        if verification:
            vertices = arc[(arc >= lower) & (arc <= upper)]
            samples = np.unique(np.r_[samples, vertices, np.arange(lower, upper, .025)])
        slopes.append(cosines(points, arc, samples))
    return np.concatenate(slopes), spacing


def intersection_kind(shape):
    if shape.is_empty:
        return ()
    if hasattr(shape, "geoms"):
        return tuple(sorted(kind for part in shape.geoms for kind in intersection_kind(part)))
    return (shape.geom_type,)


class GeometricConstraint:
    def __init__(self, pipeline, stage):
        self.pipeline = pipeline
        self.stage = stage
        self.original = sharp_geometry()
        self.editable = editable_geometry(pipeline.grid, pipeline.free)
        segments = []
        for polygon in getattr(self.editable, "geoms", (self.editable,)):
            for ring in (polygon.exterior, *polygon.interiors):
                points = list(ring.coords)
                segments.extend(LineString(pair) for pair in zip(points, points[1:]))
        self.clip_segments = segments

    def sample(self, design):
        s = self.stage
        density = self.pipeline.evaluate(design, s.q, s.beta, s.alpha_max)
        raw = reconstruct(self.pipeline.grid, density.projected)
        clipped = _clipped(self.original, self.editable, raw)
        chains = wall_chains(clipped)
        cosine, spacing = sample_angles(chains)
        # Detect marching-squares sign changes, ambiguous saddle switches, and
        # contour crossings that move between segments of the fixed clipping mask.
        field = density.projected.reshape(self.pipeline.grid.shape) - .5
        a, b, c, d = field[:-1, :-1], field[:-1, 1:], field[1:, 1:], field[1:, :-1]
        ambig = ((a >= 0) == (c >= 0)) & ((b >= 0) == (d >= 0)) & ((a >= 0) != (b >= 0))
        saddle = np.sign((a * c - b * d)[ambig]).tobytes()
        contour = MultiLineString(raw.contour_segments_mm)
        clipping = tuple(intersection_kind(contour.intersection(edge)) for edge in self.clip_segments)
        signature = ((field >= 0).tobytes(), saddle, clipping, raw.port_components,
                     tuple(len(p.interiors) for p in raw.polygons), clipped.port_components)
        return float(TAU * logsumexp((LIMIT - cosine) / TAU) / SCALE), signature, spacing, clipped, chains

    def value(self, design, signature=None):
        sampled = self.sample(design)
        if signature is not None and sampled[1] != signature:
            raise UnresolvedGeometry("Density perturbation changed threshold-contour or clipping signature")
        return sampled[0]

    def gradient(self, design, base, step):
        value, signature = base[:2]
        result = np.empty(len(design))
        trial = design.copy()
        for index, x in enumerate(design):
            try:
                roundoff = 64 * np.finfo(float).eps
                if roundoff <= x <= 1. - roundoff:
                    span = min(step, x, 1. - x)
                    trial[index] = x + span
                    plus = self.value(trial, signature)
                    trial[index] = x - span
                    result[index] = (plus - self.value(trial, signature)) / (2. * span)
                else:
                    signed = step if x < .5 else -step
                    trial[index] = x + signed
                    result[index] = (self.value(trial, signature) - value) / signed
            except UnresolvedGeometry as error:
                raise UnresolvedGeometry(f"Control {index}, density {x:.12g}, step {step:g}: {error}") from error
            finally:
                trial[index] = x
        return result
