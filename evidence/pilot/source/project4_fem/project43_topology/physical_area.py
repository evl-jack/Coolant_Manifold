"""Differentiable upper bound on the area reconstructed at density 0.5.

The fixed collars contribute exactly 41 mm². Their interiors and the design
rectangle are disjoint. Inside the design rectangle we use the same padded
cell-center patches as marching squares, including the half-width edge patches.

For each patch let f_i be density minus 0.5. At an alternating-sign patch,
the fluid-connected marching-squares polygon contains both the disconnected
branch and the exact-saddle polygons. The disconnected fluid corner triangles
avoid the connected branch's two solid corner triangles. At an exact saddle,
the edge intersections share the saddle's two axis coordinates, and both
fluid-to-saddle rings avoid the same solid triangles. Use this upper branch
for every alternating pattern to define an upper reconstruction.

The upper reconstruction is monotone in each corner value. For fixed signs,
increasing a corner moves each affected crossing toward the solid side. The
affected corner triangles expand or shrink accordingly. For adjacent signs,
the two endpoints of the separating chord move on opposite edges and the
fluid-side span is their linear interpolation. When a corner changes from
solid to fluid, the old and new polygons agree at the limiting endpoint on
that corner, except when a second, opposite fluid corner appears. In that
case the upper branch adds the central region and preserves the old fluid
triangle. A zero corner is treated as fluid, matching marching squares.

Let q_i be a C1 probability of retaining negative corner i as solid, equal to
one for f_i <= -transition and zero for f_i >= 0. Sum over solid subsets S
with weights prod(i in S) q_i * prod(i outside S) (1-q_i). These nonnegative
weights sum to one. Every subset with positive weight keeps its truly solid
corners unchanged and lifts the other corners to a smooth upper positive
part. For an edge from solid magnitude n to fluid magnitude p, replace the
solid reach n/(n+p) by n/(n+p+regularization), a further enlargement of the
fluid region. Each term thus bounds the same original polygon from above.
The derivative differentiates the subset weight and the crossing reaches,
then accumulates repeated boundary samples back to the design vector. At
f_i = 0 both the weight and its derivative vanish for every term retaining
corner i as solid, so the clipped negative magnitude introduces no kink.
"""

from itertools import product

import numpy as np

from .density import FIXED_FLUID_AREA_MM2


class ReconstructedAreaUpperBound:
    """Smooth upper bound with analytic gradient in projected density."""

    def __init__(self, grid, *, transition=.02, regularization=1e-6):
        if not np.isfinite(transition) or not 0. < transition <= .5:
            raise ValueError("The solid transition must be in (0, 0.5]")
        if not np.isfinite(regularization) or regularization <= 0.:
            raise ValueError("Crossing regularization must be positive and finite")
        self.size = grid.size
        self.transition = float(transition)
        self.regularization = float(regularization)

        centers = grid.centers_mm.reshape(*grid.shape, 2)
        x = np.r_[grid.x_edges_mm[0], centers[0, :, 0], grid.x_edges_mm[-1]]
        y = np.r_[grid.y_edges_mm[0], centers[:, 0, 1], grid.y_edges_mm[-1]]
        self.patch_area = np.outer(np.diff(y), np.diff(x)).ravel()
        indices = np.pad(np.arange(grid.size).reshape(grid.shape), 1, mode="edge")
        self.corner_indices = np.stack((indices[:-1, :-1], indices[:-1, 1:],
                                        indices[1:, 1:], indices[1:, :-1]), axis=-1).reshape(-1, 4)

    def evaluate(self, projected):
        """Return the upper area in mm² and its projected-density gradient."""
        projected = np.asarray(projected, dtype=float)
        if projected.shape != (self.size,) or not np.isfinite(projected).all():
            raise ValueError("Expected one finite projected density value per design cell")
        offset = projected[self.corner_indices] - .5
        t = np.clip(-offset / self.transition, 0., 1.)
        solid_probability = t**2 * (3. - 2. * t)
        probability_derivative = -6. * t * (1. - t) / self.transition

        norm = np.hypot(offset, self.regularization)
        positive = .5 * (offset + norm) + self.regularization
        positive_derivative = .5 * (1. + offset / norm)
        negative = np.maximum(-offset, 0.)
        negative_derivative = -(offset < 0.).astype(float)

        # Entry (i,j) is the solid-side crossing reach for a solid corner i
        # and a fluid corner j. The two derivatives use the same indexing.
        n, p = negative[:, :, None], positive[:, None, :]
        denominator = n + p
        reach = n / denominator
        solid_derivative = (p / denominator) * (negative_derivative[:, :, None] / denominator)
        fluid_derivative = -reach * (positive_derivative[:, None, :] / denominator)

        fraction = np.zeros(len(offset))
        gradient = np.zeros_like(offset)
        for pattern in product((False, True), repeat=4):
            solid = np.asarray(pattern)
            factors = np.where(solid, solid_probability, 1. - solid_probability)
            weight = factors.prod(axis=1)
            active = weight > 0.
            if not active.any() or solid.all():
                continue
            local_factors = factors[active]
            weight_derivative = np.column_stack([
                local_factors[:, [j for j in range(4) if j != i]].prod(axis=1)
                * probability_derivative[active, i] * (1. if solid[i] else -1.)
                for i in range(4)
            ])
            area, area_derivative = self._patch_fraction(
                solid, reach[active], solid_derivative[active], fluid_derivative[active])
            fraction[active] += weight[active] * area
            gradient[active] += weight_derivative * area[:, None] + weight[active, None] * area_derivative

        gradient *= self.patch_area[:, None]
        return (float(FIXED_FLUID_AREA_MM2 + self.patch_area @ fraction),
                np.bincount(self.corner_indices.ravel(), weights=gradient.ravel(), minlength=self.size))

    @staticmethod
    def _patch_fraction(solid, reach, solid_derivative, fluid_derivative):
        derivative = np.zeros((len(reach), 4))
        solid_corners = np.flatnonzero(solid)
        if len(solid_corners) == 0:
            return np.ones(len(reach)), derivative

        if len(solid_corners) == 3:
            f = int(np.flatnonzero(~solid)[0])
            a, b = (f - 1) % 4, (f + 1) % 4
            # A single fluid triangle uses the complement of each solid reach.
            # Compute the equivalent derivative by differentiating its two legs.
            leg_a = 1. - reach[:, a, f]
            leg_b = 1. - reach[:, b, f]
            derivative[:, a] = -.5 * leg_b * solid_derivative[:, a, f]
            derivative[:, b] = -.5 * leg_a * solid_derivative[:, b, f]
            derivative[:, f] = -.5 * (leg_b * fluid_derivative[:, a, f]
                                      + leg_a * fluid_derivative[:, b, f])
            return .5 * leg_a * leg_b, derivative

        if len(solid_corners) == 2 and np.diff(solid_corners)[0] != 2:
            fraction = np.ones(len(reach))
            for s in solid_corners:
                a, b = (s - 1) % 4, (s + 1) % 4
                f = b if solid[a] else a
                fraction -= .5 * reach[:, s, f]
                derivative[:, s] -= .5 * solid_derivative[:, s, f]
                derivative[:, f] -= .5 * fluid_derivative[:, s, f]
            return fraction, derivative

        # A lone solid corner and the upper alternating branch both remove
        # one solid triangle per solid corner.
        fraction = np.ones(len(reach))
        for s in solid_corners:
            a, b = (s - 1) % 4, (s + 1) % 4
            leg_a, leg_b = reach[:, s, a], reach[:, s, b]
            fraction -= .5 * leg_a * leg_b
            derivative[:, s] -= .5 * (leg_b * solid_derivative[:, s, a]
                                      + leg_a * solid_derivative[:, s, b])
            derivative[:, a] -= .5 * leg_b * fluid_derivative[:, s, a]
            derivative[:, b] -= .5 * leg_a * fluid_derivative[:, s, b]
        return fraction, derivative
