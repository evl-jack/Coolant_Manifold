"""Density-field tools for Brinkman-RANS experiments on a Cartesian grid."""

import numpy as np


def resistance(gamma, alpha_max, power=2):
    """Return alpha_max*(1-gamma)**power for gamma=1 fluid, 0 solid."""
    if alpha_max <= 0 or power <= 0:
        raise ValueError("alpha_max and power must be positive")
    return alpha_max * np.clip(1 - np.asarray(gamma), 0, 1) ** power


def momentum_source(gamma, u, v, alpha_max, power=2):
    """Return the per-volume Brinkman sources for the two momentum equations."""
    alpha = resistance(gamma, alpha_max, power)
    return -alpha * np.asarray(u), -alpha * np.asarray(v)


def _edge_crossing(p0, p1, g0, g1, threshold, method, extra=None):
    if (g0 - threshold) * (g1 - threshold) >= 0:
        return None
    linear = p0 + (threshold - g0) * (p1 - p0) / (g1 - g0)
    if method == "linear" or extra is None:
        return linear
    px, gx = extra
    values = np.asarray(gx)
    if np.any(np.isclose(np.diff(values), 0)) or np.any(np.diff(values) * np.diff(values)[0] < 0):
        return linear
    coeff = np.polynomial.polynomial.polyfit(px, values, 2)
    roots = np.polynomial.polynomial.polyroots([coeff[0] - threshold, coeff[1], coeff[2]])
    roots = roots[np.abs(roots.imag) < 1e-10].real
    roots = roots[(roots >= min(p0, p1)) & (roots <= max(p0, p1))]
    return roots[np.argmin(np.abs(roots - linear))] if len(roots) else linear


def extract_wall_segments(x, y, gamma, threshold=.5, method="linear"):
    """Extract interface segments from cell-center density on a Cartesian grid.

    Returns an array shaped (n_segments, 2, 2), in the same coordinate units as
    x and y. Quadratic interpolation falls back to linear at nonmonotone or
    invalid crossings. Four-crossing saddle cells use a deterministic pairing.
    """
    x, y, gamma = np.asarray(x), np.asarray(y), np.asarray(gamma)
    if gamma.shape != (len(x), len(y)):
        raise ValueError("gamma must have shape (len(x), len(y))")
    if method not in ("linear", "quadratic"):
        raise ValueError("method must be 'linear' or 'quadratic'")
    segments = []

    def cross(i, j, axis):
        if axis == 0:
            p0, p1, g0, g1 = x[i], x[i + 1], gamma[i, j], gamma[i + 1, j]
            if method == "quadratic" and len(x) >= 3:
                ids = (i - 1, i, i + 1) if i > 0 else (i, i + 1, i + 2)
                extra = (x[list(ids)], gamma[list(ids), j])
            else:
                extra = None
            px = _edge_crossing(p0, p1, g0, g1, threshold, method, extra)
            return None if px is None else np.array([px, y[j]])
        p0, p1, g0, g1 = y[j], y[j + 1], gamma[i, j], gamma[i, j + 1]
        if method == "quadratic" and len(y) >= 3:
            ids = (j - 1, j, j + 1) if j > 0 else (j, j + 1, j + 2)
            extra = (y[list(ids)], gamma[i, list(ids)])
        else:
            extra = None
        py = _edge_crossing(p0, p1, g0, g1, threshold, method, extra)
        return None if py is None else np.array([x[i], py])

    for i in range(len(x) - 1):
        for j in range(len(y) - 1):
            edges = [cross(i, j, 0), cross(i + 1, j, 1),
                     cross(i, j + 1, 0), cross(i, j, 1)]
            hits = [p for p in edges if p is not None]
            if len(hits) == 2:
                segments.append(hits)
            elif len(hits) == 4:
                segments.extend(((edges[0], edges[1]), (edges[2], edges[3])))
    return np.asarray(segments, dtype=float).reshape(-1, 2, 2)


def wall_distance(points, segments):
    """Compute each point's shortest Euclidean distance to line segments."""
    points, segments = np.asarray(points, dtype=float), np.asarray(segments, dtype=float)
    if not len(segments):
        raise ValueError("wall segments are empty")
    starts, vectors = segments[:, 0], segments[:, 1] - segments[:, 0]
    t = np.sum((points[:, None] - starts) * vectors, axis=2) / np.maximum(np.sum(vectors ** 2, axis=1), 1e-30)
    nearest = starts + np.clip(t, 0, 1)[:, :, None] * vectors
    return np.min(np.linalg.norm(points[:, None] - nearest, axis=2), axis=1)
