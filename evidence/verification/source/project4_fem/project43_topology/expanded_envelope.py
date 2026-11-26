"""Compare the archived design box with a downward-only full-density expansion."""

import argparse
import os
from pathlib import Path

import numpy as np

from ..fem import TaylorHood
from ..junction import boundary_data, strip_triangles
from .density import DesignGrid, DensityPipeline, FIXED_FLUID_AREA_MM2
from .domain import REFERENCE_HEIGHT_MM, TopologyDomain, _refined_axis, build_domain
from .evaluator import DesignTransfer, baseline_reference
from .gate2 import _baseline_record, _continuation, environment, write_report
from .gate3_screen import fitted_reference
from .reconstruction import _PORTS_MM


ADDED_DEPTH_MM = 1.0
ORIGINAL_FLUID_AREA_MM2 = 83.0


class ExpandedDesignGrid(DesignGrid):
    """Prepend design rows below y=0 without altering the original T or its collars."""

    def __init__(self, spacing_mm):
        super().__init__(spacing_mm)
        rows = ADDED_DEPTH_MM / self.spacing_mm
        count = round(rows)
        if count < 1 or not np.isclose(rows, count, rtol=0.0, atol=1e-10):
            raise ValueError("The design spacing must align exactly with the added strip")
        below_wall = np.linspace(-ADDED_DEPTH_MM, 0.0, count + 1)[:-1]
        self.original_shape = self.shape
        self.added_rows = count
        self.y_edges_mm = np.r_[below_wall, self.y_edges_mm]
        self.shape = (self.shape[0] + count, self.shape[1])
        self.size = self.shape[0] * self.shape[1]
        self.base_raw = np.r_[np.zeros(count * self.shape[1]), self.base_raw]
        centers_y = (self.y_edges_mm[:-1] + self.y_edges_mm[1:]) / 2
        centers_x = (self.x_edges_mm[:-1] + self.x_edges_mm[1:]) / 2
        xx, yy = np.meshgrid(centers_x, centers_y, indexing="xy")
        self.centers_mm = np.column_stack((xx.ravel(), yy.ravel()))


def _check_domain(grid, domain, original_fluid):
    fem = domain.fem
    area = fem.weights.sum(axis=1) * REFERENCE_HEIGHT_MM**2
    inside = domain.cell_index >= 0
    by_cell = np.bincount(domain.cell_index[inside], weights=area[inside], minlength=grid.size)
    np.testing.assert_allclose(area.sum(), FIXED_FLUID_AREA_MM2 + grid.size * grid.cell_area_mm2, rtol=1e-12)
    np.testing.assert_allclose(area[~inside].sum(), FIXED_FLUID_AREA_MM2, rtol=1e-12)
    np.testing.assert_allclose(area[original_fluid].sum(), ORIGINAL_FLUID_AREA_MM2, rtol=1e-12)
    np.testing.assert_allclose(by_cell, grid.cell_area_mm2, rtol=1e-12)
    np.testing.assert_array_equal(grid.base_raw[domain.cell_index[inside]], original_fluid[inside])

    iy, ix = np.unravel_index(domain.cell_index[inside], grid.shape)
    corners = fem.points[fem.triangles[inside]] * REFERENCE_HEIGHT_MM
    lower = np.column_stack((grid.x_edges_mm[ix], grid.y_edges_mm[iy]))
    upper = np.column_stack((grid.x_edges_mm[ix + 1], grid.y_edges_mm[iy + 1]))
    if (np.any(corners.min(axis=1) < lower - 1e-12)
            or np.any(corners.max(axis=1) > upper + 1e-12)):
        raise ValueError("An expanded-envelope triangle crosses a design-cell face")

    for edges, (start, stop) in zip(domain.ports, _PORTS_MM):
        points = fem.points[edges[:, :2]] * REFERENCE_HEIGHT_MM
        np.testing.assert_allclose(points.min(axis=(0, 1)), start, rtol=0.0, atol=1e-12)
        np.testing.assert_allclose(points.max(axis=(0, 1)), stop, rtol=0.0, atol=1e-12)
        width = np.linalg.norm(points[:, 1] - points[:, 0], axis=1).sum()
        np.testing.assert_allclose(width, np.linalg.norm(np.subtract(stop, start)), rtol=1e-12)


def build_expanded_domain(grid, max_spacing_mm):
    """Use the verified FEM on the original T union the expanded design box."""
    if not isinstance(grid, ExpandedDesignGrid):
        raise TypeError("Use ExpandedDesignGrid with the expanded-envelope domain")
    if not np.isfinite(max_spacing_mm) or max_spacing_mm <= 0.0:
        raise ValueError("The maximum FEM axis spacing must be finite and positive")
    x_axis = _refined_axis(np.r_[0.0, grid.x_edges_mm, 9.5, 10.5, 20.0], max_spacing_mm)
    y_axis = _refined_axis(np.r_[grid.y_edges_mm, 0.0, 4.0, 7.0], max_spacing_mm)
    xx, yy = np.meshgrid(x_axis, y_axis)
    points_mm = np.column_stack((xx.ravel(), yy.ravel()))
    triangles = strip_triangles(np.arange(xx.size).reshape(xx.shape))
    x, y = points_mm[triangles].mean(axis=1).T
    original_fluid = ((0.0 < y) & (y < 4.0)) | ((9.5 < x) & (x < 10.5) & (4.0 < y))
    designable = ((grid.x_edges_mm[0] < x) & (x < grid.x_edges_mm[-1])
                  & (grid.y_edges_mm[0] < y) & (y < grid.y_edges_mm[-1]))
    keep = original_fluid | designable
    used_vertices, local_vertices = np.unique(triangles[keep], return_inverse=True)
    fem = TaylorHood(points_mm[used_vertices] / REFERENCE_HEIGHT_MM, local_vertices.reshape(-1, 3))

    cell_index = np.full(len(fem.triangles), -1, dtype=np.intp)
    ix = np.searchsorted(grid.x_edges_mm, x[designable], side="right") - 1
    iy = np.searchsorted(grid.y_edges_mm, y[designable], side="right") - 1
    cell_index[designable[keep]] = np.ravel_multi_index((iy, ix), grid.shape)
    fixed, prescribed, ports, flux = boundary_data(fem)
    domain = TopologyDomain(fem, fixed, prescribed, ports, flux, cell_index, max_spacing_mm)
    _check_domain(grid, domain, original_fluid[keep])
    return domain


def _baseline(grid, domain, alpha_max, filter_radius_mm, betas):
    pipeline = DensityPipeline(grid, filter_radius_mm=filter_radius_mm)
    raw_area = FIXED_FLUID_AREA_MM2 + grid.cell_area_mm2 * grid.base_raw.sum()
    np.testing.assert_allclose(raw_area, ORIGINAL_FLUID_AREA_MM2, rtol=1e-12)
    expanded_count = getattr(grid, "added_rows", 0) * grid.shape[1]
    density_states = {}
    density_records = {}
    for beta in betas:
        state = pipeline.evaluate(grid.base_raw, q=0.01, beta=beta, alpha_max=alpha_max)
        density_states[beta] = state
        density_records[str(beta)] = {
            "filtered_total_mm2": FIXED_FLUID_AREA_MM2 + grid.cell_area_mm2 * state.filtered.sum(),
            "projected_total_mm2": state.area_mm2,
            "projected_margin_to_83_mm2": ORIGINAL_FLUID_AREA_MM2 - state.area_mm2,
            "filtered_in_added_strip_mm2": grid.cell_area_mm2 * state.filtered[:expanded_count].sum(),
            "projected_in_added_strip_mm2": grid.cell_area_mm2 * state.projected[:expanded_count].sum(),
        }
    np.testing.assert_allclose(density_records["0.0"]["filtered_total_mm2"], raw_area, rtol=1e-12)

    # At beta=0, solve both the binary sharp-wall penalty mask and the setup's
    # usual filtered initialization to expose the two distinct envelope effects.
    transfer = DesignTransfer(domain.cell_index, grid.size)
    binary = _continuation(domain.fem, domain.fixed, domain.prescribed,
                           transfer.forward(alpha_max * (1.0 - grid.base_raw)))
    _, filtered = baseline_reference(domain, pipeline, alpha_max)
    report = {
        "bounds_mm": [[grid.x_edges_mm[0], grid.x_edges_mm[-1]],
                      [grid.y_edges_mm[0], grid.y_edges_mm[-1]]],
        "shape_yx": grid.shape,
        "design_cells": grid.size,
        "design_area_mm2": grid.size * grid.cell_area_mm2,
        "raw_fluid_area_mm2": raw_area,
        "fixed_fluid_area_mm2": FIXED_FLUID_AREA_MM2,
        "computational_area_mm2": FIXED_FLUID_AREA_MM2 + grid.size * grid.cell_area_mm2,
        "density_budget": density_records,
        "unfiltered_binary_penalty": _baseline_record(domain, binary),
        "filtered_beta_zero_penalty": _baseline_record(domain, filtered),
    }
    return report, density_states


def matched_baseline(design_spacing_mm, analysis_spacing_mm, alpha_max, filter_radius_mm, betas):
    """Return separate same-geometry baseline solves before comparing optimized designs."""
    sol_job_id = os.environ.get("SLURM_JOB_ID")
    if not sol_job_id:
        raise RuntimeError("Run this numerical diagnostic in a Sol Slurm job")
    betas = tuple(sorted({0.0, *(float(beta) for beta in betas)}))
    if not np.all(np.isfinite(betas)) or min(betas) < 0.0:
        raise ValueError("Projection betas must be finite and nonnegative")
    original = DesignGrid(design_spacing_mm)
    expanded = ExpandedDesignGrid(design_spacing_mm)
    np.testing.assert_array_equal(expanded.base_raw.reshape(expanded.shape)[expanded.added_rows:],
                                  original.base_raw.reshape(original.shape))
    original_report, original_states = _baseline(original, build_domain(original, analysis_spacing_mm),
                                                 alpha_max, filter_radius_mm, betas)
    expanded_report, expanded_states = _baseline(expanded, build_expanded_domain(expanded, analysis_spacing_mm),
                                                  alpha_max, filter_radius_mm, betas)
    fitted = fitted_reference("sharp", analysis_spacing_mm)
    deltas = {}
    for beta in betas:
        restricted = expanded_states[beta].projected.reshape(expanded.shape)[expanded.added_rows:].ravel()
        difference = restricted - original_states[beta].projected
        deltas[str(beta)] = {
            "common_cells_projected_density_rms": np.sqrt(np.mean(difference**2)),
            "common_cells_projected_density_max_abs": np.max(np.abs(difference)),
            "common_cells_projected_fluid_volume_change_mm2": original.cell_area_mm2 * difference.sum(),
            "total_projected_fluid_volume_change_mm2": expanded_states[beta].area_mm2
                                                       - original_states[beta].area_mm2,
        }
    for envelope in (original_report, expanded_report):
        for mode in ("unfiltered_binary_penalty", "filtered_beta_zero_penalty"):
            row = envelope[mode]
            row["split_difference_from_fitted_pp"] = 100 * (row["side_fraction"] - fitted["side_fraction"])
            row["power_relative_difference_from_fitted"] = row["power_w_per_m"] / fitted["power_w_per_m"] - 1.0
    return {
        "sol_job_id": sol_job_id,
        "environment": environment(),
        "design_spacing_mm": design_spacing_mm,
        "analysis_spacing_mm": analysis_spacing_mm,
        "alpha_max": alpha_max,
        "filter_radius_mm": filter_radius_mm,
        "same_original_sharp_fitted_reference": fitted,
        "original": original_report,
        "expanded": expanded_report,
        "expanded_minus_original": deltas,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="Fresh output directory")
    parser.add_argument("--design-spacing-mm", type=float, default=0.5)
    parser.add_argument("--analysis-spacing-mm", type=float, default=0.5)
    parser.add_argument("--alpha-max", type=float, default=1e5)
    parser.add_argument("--filter-radius-mm", type=float, default=0.5)
    parser.add_argument("--projection-betas", nargs="+", type=float, default=[2.0, 8.0])
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical diagnostic in a Sol Slurm job")
    args.output.mkdir(parents=True)
    if not write_report(args.output, "matched_envelope_baseline", lambda: matched_baseline(
        args.design_spacing_mm, args.analysis_spacing_mm, args.alpha_max, args.filter_radius_mm, args.projection_betas,
    )):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
