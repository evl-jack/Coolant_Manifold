"""Render matched PNG evidence from archived fixed-pipe pilot states."""

from pathlib import Path

import numpy as np

from .domain import REFERENCE_HEIGHT_MM
from .evaluator import HEIGHT_M, VELOCITY_SCALE, VISCOSITY


FULL_VIEW = (-.45, 20.45, -.38, 7.38)
JUNCTION_VIEW = (7.35, 12.65, 3.3, 7.25)
FLUID = "#216181"
ORIGINAL = "#aa5b20"
OPTIMIZED = "#147d83"
FIXED = "#dde7eb"


def plot_stage(output, domain, pipeline, baseline, candidate, history):
    """Save the actual immersed fields and threshold geometry, without a CFD replay."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import PathPatch, Rectangle
    from matplotlib.path import Path as PlotPath
    from matplotlib.tri import Triangulation
    from shapely.geometry import box
    from shapely.geometry.polygon import orient

    target = Path(output)
    target.mkdir(parents=True, exist_ok=True)
    names = []
    grid = pipeline.grid
    free = np.asarray(pipeline.free, dtype=bool).reshape(grid.shape)
    rows, columns = np.nonzero(free)
    xmin, xmax = grid.x_edges_mm[[columns.min(), columns.max() + 1]]
    ymin, ymax = grid.y_edges_mm[[rows.min(), rows.max() + 1]]
    edit_box = box(xmin, ymin, xmax, ymax)
    original = baseline["geometry"].geometry
    optimized = candidate["geometry"].geometry
    fixed_fluid = original.difference(edit_box)

    def save(fig, name):
        filename = f"{name}.png"
        try:
            fig.savefig(target / filename, dpi=175, facecolor="white")
        finally:
            plt.close(fig)
        names.append(filename)

    def parts(geometry):
        if geometry.is_empty:
            return
        if geometry.geom_type == "Polygon":
            yield geometry
        elif hasattr(geometry, "geoms"):
            for part in geometry.geoms:
                yield from parts(part)

    def fill(ax, geometry, color):
        for polygon in parts(geometry):
            polygon = orient(polygon, sign=1.)
            rings = [PlotPath(np.asarray(ring.coords), closed=True)
                     for ring in (polygon.exterior, *polygon.interiors)]
            ax.add_patch(PathPatch(PlotPath.make_compound_path(*rings), facecolor=color, edgecolor="none"))

    def trace(ax, geometry, color, *, dashed=False, width=1.1):
        for polygon in parts(geometry):
            for ring in (polygon.exterior, *polygon.interiors):
                ax.plot(*ring.xy, color=color, linewidth=width, linestyle="--" if dashed else "-")

    def geometry_axis(ax, title, view=FULL_VIEW):
        ax.set(title=title, aspect="equal", xlim=view[:2], ylim=view[2:], xlabel="x (mm)", ylabel="y (mm)")
        ax.tick_params(labelsize=8, color="#79848c")
        for spine in ax.spines.values():
            spine.set_color("#a4aeb4")

    def rectangle(ax, *, face="none", edge="#a44f24", dashed=True):
        ax.add_patch(Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, facecolor=face,
                               edgecolor=edge, linewidth=1.1, linestyle="--" if dashed else "-"))

    def render_overlay(ax):
        fill(ax, original.intersection(optimized), FIXED)
        fill(ax, optimized.difference(original), "#a9d8d6")
        fill(ax, original.difference(optimized), "#f4d1a5")
        trace(ax, optimized, OPTIMIZED, width=1.35)
        trace(ax, original, ORIGINAL, dashed=True)
        ax.legend(handles=[Line2D([0], [0], color=ORIGINAL, linestyle="--", label="Original wall"),
                           Line2D([0], [0], color=OPTIMIZED, label="Reconstructed wall")],
                  loc="lower right", fontsize=7, frameon=True)

    def render_density(ax, field):
        fill(ax, fixed_fluid, FIXED)
        trace(ax, original, "#72858e", dashed=True, width=.65)
        data = np.ma.masked_where(~free, np.asarray(field).reshape(grid.shape))
        image = ax.pcolormesh(grid.x_edges_mm, grid.y_edges_mm, data, shading="flat",
                              cmap="viridis", vmin=0., vmax=1., rasterized=True)
        rectangle(ax, edge="#61757f", dashed=False)
        return image

    for name, title, geometry in (
        ("original_physical_geometry", "Original physical geometry", original),
        ("optimized_physical_geometry", "Optimized physical geometry at density 0.5", optimized),
    ):
        fig, ax = plt.subplots(figsize=(9.3, 4.1), layout="constrained")
        geometry_axis(ax, title)
        fill(ax, geometry, FLUID)
        trace(ax, geometry, "#154b67", width=.8)
        save(fig, name)

    fig, ax = plt.subplots(figsize=(9.3, 4.1), layout="constrained")
    geometry_axis(ax, f"Brinkman transition region · {free.sum()} editable cells")
    fill(ax, fixed_fluid, FIXED)
    rectangle(ax, face="#ffe4c5", edge=ORIGINAL)
    trace(ax, original, "#526b7a", dashed=True)
    save(fig, "editable_brinkman_region")

    for name, title, field in (
        ("initial_density", "Initial design density", baseline["density"].raw),
        ("optimized_density", "Optimized projected density", candidate["density"].projected),
    ):
        fig, ax = plt.subplots(figsize=(9.3, 4.1), layout="constrained")
        geometry_axis(ax, title)
        fig.colorbar(render_density(ax, field), ax=ax, shrink=.75, label="Fluid density")
        save(fig, name)

    fig, ax = plt.subplots(figsize=(9.3, 4.1), layout="constrained")
    geometry_axis(ax, "Before and after physical geometry")
    render_overlay(ax)
    save(fig, "geometry_overlay")

    fem = domain.fem
    a, b, c, ab, bc, ca = fem.velocity_dofs.T
    subtriangles = np.vstack([np.column_stack(nodes) for nodes in (
        (a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca),
    )])
    velocities = Triangulation(*(fem.velocity_points * REFERENCE_HEIGHT_MM).T, triangles=subtriangles)
    pressure = Triangulation(*(fem.points * REFERENCE_HEIGHT_MM).T, triangles=fem.triangles)

    def speeds(state):
        return VELOCITY_SCALE * np.linalg.norm(state[:2 * fem.n_u].reshape(2, fem.n_u), axis=0)

    def pressures(state):
        return VISCOSITY * VELOCITY_SCALE / HEIGHT_M * state[2 * fem.n_u:]

    states = [np.asarray(case["state"]) for case in (baseline, candidate)]
    for name, label, mesh, field, cmap in (
        ("velocity_magnitude_before_after", "Velocity magnitude (m/s)", velocities, speeds, "magma"),
        ("pressure_before_after", "Pressure (Pa)", pressure, pressures, "cividis"),
    ):
        values = [field(state) for state in states]
        low = 0. if field is speeds else min(float(data.min()) for data in values)
        high = max(float(data.max()) for data in values)
        if low == high:
            high = low + max(abs(low), 1.) * 1e-8
        levels = np.linspace(low, high, 61)
        fig, axes = plt.subplots(2, 1, figsize=(9.4, 7.), layout="constrained")
        fig.suptitle(f"Immersed Navier–Stokes · {label}", fontsize=13)
        for ax, geometry, title, data in zip(axes, (original, optimized), ("Before", "After"), values):
            geometry_axis(ax, title)
            image = ax.tricontourf(mesh, data, levels=levels, cmap=cmap)
            trace(ax, geometry, "#dae1e4", dashed=True, width=.65)
        fig.colorbar(image, ax=axes, shrink=.8, label=label)
        save(fig, name)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6), layout="constrained")
    geometry_axis(axes[0], "Original and reconstructed walls", JUNCTION_VIEW)
    render_overlay(axes[0])
    geometry_axis(axes[1], "Projected density in the transition", JUNCTION_VIEW)
    fig.colorbar(render_density(axes[1], candidate["density"].projected), ax=axes[1], shrink=.75,
                 label="Fluid density")
    save(fig, "junction_magnified")

    seen = set()
    evaluations = []
    for entry in history:
        key = (entry["stage"], entry["phase"], entry["unique_evaluations"])
        if key not in seen:
            evaluations.append(entry)
            seen.add(key)
    for name, label, field, scale, reference in (
        ("pumping_power_history", "Normalized pumping power", "normalized_power", 1., None),
        ("side_flow_fraction_history", "Side-outlet flow (%)", "side_fraction", 100., "reference_side_fraction"),
    ):
        fig, ax = plt.subplots(figsize=(9.1, 3.5), layout="constrained")
        x = np.arange(1, len(evaluations) + 1)
        ax.plot(x, [scale * item[field] for item in evaluations], "o-", color=FLUID, markersize=2., linewidth=1.15,
                label="Evaluated design")
        if reference:
            ax.step(x, [scale * item[reference] for item in evaluations], where="post", color=ORIGINAL,
                    linestyle="--", linewidth=1.15, label="Original design at each stage")
        for index, entry in enumerate(evaluations):
            if index == 0 or entry["stage"] != evaluations[index - 1]["stage"]:
                ax.axvline(x[index], color="#b3bdc2", linewidth=.7, zorder=0)
                ax.annotate(f"Stage {entry['stage']}", (x[index], .98), xycoords=("data", "axes fraction"),
                            xytext=(4, -2), textcoords="offset points", fontsize=7.5, va="top", color="#4b626e")
        ax.set(xlabel="Unique design evaluations by phase", ylabel=label)
        ax.grid(alpha=.18)
        if reference:
            ax.legend(loc="best", fontsize=8)
        save(fig, name)
    return names
