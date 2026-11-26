"""Generate the fixed-pipe Stokes versus Navier–Stokes evidence and figures."""

import argparse
import json
from pathlib import Path
from zipfile import ZipFile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as PlotPath
from matplotlib.tri import Triangulation
import numpy as np
from shapely.geometry import shape
from shapely.geometry.polygon import orient


DESIGNS = ("original", "stokes_optimized", "navier_stokes_optimized")
DESIGN_LABELS = ("Original", "Stokes optimized", "Navier–Stokes optimized")
MODELS = ("stokes", "navier_stokes")
MODEL_LABELS = ("Stokes", "Navier–Stokes")
COLORS = ("#8b969c", "#9d512c", "#177383")
FULL_VIEW = (-.35, 20.35, -.30, 7.3)
JUNCTION_VIEW = (7.55, 12.45, 3.45, 7.18)
FIGURE_DIR = "figures"


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def unwrap(data):
    return data["results"] if "completed" in data else data


def model_data(pilot, kkt_path):
    config = unwrap(read_json(pilot / "configuration.json"))
    stages = [unwrap(read_json(path)) for path in sorted(pilot.glob("stage_*/summary.json"))]
    with (pilot / "progress.jsonl").open(encoding="utf-8") as stream:
        progress = [json.loads(line) for line in stream if line.strip()]
    if kkt_path is not None and kkt_path.suffix == ".zip":
        with ZipFile(kkt_path) as report:
            kkt_stages = unwrap(json.loads(report.read("kkt-results/fixed_pipe_kkt.json")))["stages"]
    elif kkt_path is not None:
        audit = unwrap(read_json(kkt_path))
        kkt_stages = audit["stages"] if "stages" in audit else [audit]
    else:
        kkt_stages = [unwrap(read_json(path)) for path in sorted(pilot.glob("stage_*/kkt.json"))]
    return {
        "configuration": config,
        "stages": stages,
        "progress": progress,
        "kkt": kkt_stages,
        "timing": {
            field: sum(stage.get(field, 0) for stage in stages)
            for field in ("elapsed_seconds", "forward_seconds", "adjoint_seconds", "mma_unique_evaluations")
        },
    }


def axis_for_geometry(axis, title, view):
    axis.set(title=title, xlabel="x (mm)", ylabel="y (mm)", aspect="equal", xlim=view[:2], ylim=view[2:])
    axis.tick_params(labelsize=8)
    for spine in axis.spines.values():
        spine.set_color("#aab6bc")


def polygons(geometry):
    if geometry.geom_type == "Polygon":
        yield geometry
    elif hasattr(geometry, "geoms"):
        for child in geometry.geoms:
            yield from polygons(child)


def region(axis, geometry, color):
    for part in polygons(geometry):
        polygon = orient(part, sign=1.)
        paths = [PlotPath(np.asarray(ring.coords), closed=True)
                 for ring in (polygon.exterior, *polygon.interiors)]
        axis.add_patch(PathPatch(PlotPath.make_compound_path(*paths), facecolor=color, edgecolor="none"))


def contour(axis, geometry, color, **kwargs):
    for part in polygons(geometry):
        for ring in (part.exterior, *part.interiors):
            axis.plot(*ring.xy, color=color, **kwargs)


def save(figure, directory, name):
    figure.savefig(directory / f"{name}.png", dpi=215, facecolor="#ffffff")
    figure.savefig(directory / f"{name}.pdf", facecolor="#ffffff")
    plt.close(figure)


def geometry_figures(geometries, directory):
    for view, name, title in (
        (FULL_VIEW, "geometry_comparison", "Reconstructed fixed-pipe geometries"),
        (JUNCTION_VIEW, "junction_comparison", "Junction detail"),
    ):
        figure, axes = plt.subplots(1, 3, figsize=(12.2, 3.15), layout="constrained")
        figure.suptitle(title)
        for axis, geometry, label, color in zip(axes, geometries, DESIGN_LABELS, COLORS):
            axis_for_geometry(axis, label, view)
            region(axis, geometry, color)
            axis.text(.02, .04, f"Fluid area {geometry.area:.3f} mm²", transform=axis.transAxes,
                      color="#273741", fontsize=8, bbox={"facecolor": "white", "alpha": .86, "edgecolor": "none"})
        save(figure, directory, name)

    figure, axis = plt.subplots(figsize=(7.0, 4.9), layout="constrained")
    axis_for_geometry(axis, "Optimized walls at the density 0.5 contour", JUNCTION_VIEW)
    region(axis, geometries[1].intersection(geometries[2]), "#e7eef0")
    for geometry, label, color, style in zip(geometries[1:], DESIGN_LABELS[1:], COLORS[1:], ("--", "-")):
        contour(axis, geometry, color, linewidth=1.5, linestyle=style)
        axis.plot([], [], color=color, linestyle=style, label=label)
    axis.legend(loc="lower right", frameon=True, fontsize=8)
    save(figure, directory, "optimized_walls_overlay")


def density_figures(fields, geometries, directory):
    x, y = np.linspace(8, 12, 17), np.linspace(4, 6, 9)
    raw = [np.asarray(fields[(design, mode)]["design"]).reshape(8, 16)
           for design, mode in zip(DESIGNS[1:], MODELS)]
    data = (
        (raw[0], "Stokes optimized raw density", "stokes_density", "viridis", (0, 1)),
        (raw[1], "Navier–Stokes optimized raw density", "navier_stokes_density", "viridis", (0, 1)),
        (raw[0] - raw[1], "Raw density difference, Stokes minus Navier–Stokes",
         "density_difference", "RdBu_r", (-1, 1)),
    )
    for index, (values, title, filename, cmap, bounds) in enumerate(data):
        figure, axis = plt.subplots(figsize=(7., 4.6), layout="constrained")
        axis_for_geometry(axis, title, JUNCTION_VIEW)
        image = axis.pcolormesh(x, y, values, cmap=cmap, vmin=bounds[0], vmax=bounds[1], shading="flat")
        axis.add_patch(Rectangle((8, 4), 4, 2, fill=False, edgecolor="#475d69", linewidth=.7))
        if index < 2:
            contour(axis, geometries[index + 1], "#273f50", linewidth=1.)
        else:
            for geometry, color in zip(geometries[1:], COLORS[1:]):
                contour(axis, geometry, color, linewidth=.85)
        figure.colorbar(image, ax=axis, shrink=.82, label="Raw design density" if index < 2 else "Density difference")
        save(figure, directory, filename)


def plotting_field(data, key):
    if key == "pressure":
        return Triangulation(*data["pressure_points_mm"].T, triangles=data["triangles"]), data["pressure_pa"]
    a, b, c, ab, bc, ca = data["velocity_dofs"].T
    subtriangles = np.vstack([np.column_stack(nodes) for nodes in (
        (a, ab, ca), (ab, b, bc), (ca, bc, c), (ab, bc, ca),
    )])
    return (Triangulation(*data["velocity_points_mm"].T, triangles=subtriangles),
            np.linalg.norm(data["velocity_m_per_s"], axis=-1))


def field_figures(fields, geometries, directory):
    for key, label, cmap in (("velocity", "Velocity magnitude (m/s)", "magma"),
                             ("pressure", "Pressure (Pa)", "cividis")):
        plots = {case: plotting_field(field, key) for case, field in fields.items()}
        low = 0. if key == "velocity" else min(float(values.min()) for _, values in plots.values())
        high = max(float(values.max()) for _, values in plots.values())
        levels = np.linspace(low, max(high, low + max(abs(low), 1.) * 1e-9), 61)
        figure, axes = plt.subplots(2, 3, figsize=(14.0, 5.8), layout="constrained")
        for row, mode in enumerate(MODELS):
            for column, design in enumerate(DESIGNS):
                axis = axes[row, column]
                axis_for_geometry(axis, f"{DESIGN_LABELS[column]}  ·  {MODEL_LABELS[row]}", FULL_VIEW)
                grid, values = plots[(design, mode)]
                image = axis.tricontourf(grid, values, levels=levels, cmap=cmap)
                for collection in image.collections if hasattr(image, "collections") else (image,):
                    collection.set_rasterized(True)
                contour(axis, geometries[column], "#e5e7e7", linewidth=.55, linestyle="--")
        figure.colorbar(image, ax=axes, shrink=.82, label=label)
        save(figure, directory, f"{key}_comparison")


def audit_history(model):
    rows = [entry for stage in model["kkt"] for entry in stage.get("mma_history", [])]
    if rows:
        return rows
    return [dict(entry, **entry.get("kkt", {})) for entry in model["progress"]
            if "stationarity_inf" in entry or "stationarity_inf" in entry.get("kkt", {})]


def history_figures(summary, directory):
    for field, label, filename, multiplier in (
        ("normalized_power", "Normalized port pumping power", "objective_history", 1.),
        ("side_fraction", "Side-outlet flow (%)", "side_flow_history", 100.),
    ):
        figure, axes = plt.subplots(1, 2, figsize=(11., 3.7), layout="constrained", sharey=True)
        for axis, mode, label_model, color in zip(axes, MODELS, MODEL_LABELS, COLORS[1:]):
            history = [row for row in summary[mode]["progress"] if field in row]
            indices = np.arange(1, len(history) + 1)
            axis.plot(indices, [row[field] * multiplier for row in history], ".-", markersize=3., color=color)
            if field == "side_fraction":
                axis.step(indices, [row.get("reference_side_fraction", np.nan) * multiplier for row in history],
                          label="Common side-flow floor", color="#667785", linestyle="--", where="post")
                axis.legend(fontsize=7.5)
            for index, row in enumerate(history):
                if index == 0 or row.get("stage") != history[index - 1].get("stage"):
                    axis.axvline(index + 1, color="#aeb9bf", linewidth=.6, zorder=0)
                    axis.annotate(f"Stage {row.get('stage', 1)}", (index + 1, .94),
                                  xycoords=("data", "axes fraction"), xytext=(4, 0), textcoords="offset points",
                                  fontsize=7, va="top")
            axis.set(title=label_model, xlabel="Logged design evaluation", ylabel=label)
            axis.grid(alpha=.18)
        save(figure, directory, filename)

    figure, axes = plt.subplots(1, 2, figsize=(10.5, 3.7), layout="constrained")
    for axis, field, label, tol in (
        (axes[0], "stationarity_inf", "Projected stationarity residual", 1e-4),
        (axes[1], "complementarity_inf", "Complementarity residual", 1e-5),
    ):
        for mode, name, color in zip(MODELS, MODEL_LABELS, COLORS[1:]):
            rows = [entry for entry in audit_history(summary[mode]) if field in entry]
            axis.semilogy(np.arange(1, len(rows) + 1), [max(1e-16, row[field]) for row in rows],
                          color=color, linewidth=1.2, label=name)
        axis.axhline(tol, color="#667785", linestyle="--", linewidth=1., label="Tolerance")
        axis.set(xlabel="Logged audit evaluation", ylabel=label)
        axis.grid(alpha=.18)
    axes[0].legend(fontsize=7.5)
    save(figure, directory, "kkt_history")

    figure, axis = plt.subplots(figsize=(7.5, 4.2), layout="constrained")
    category = np.arange(3)
    for index, (mode, label, color) in enumerate(zip(MODELS, MODEL_LABELS, COLORS[1:])):
        value = summary[mode]["timing"]
        axis.bar(category + (index - .5) * .32, [value[k] for k in ("elapsed_seconds", "forward_seconds", "adjoint_seconds")],
                 width=.32, color=color, label=label)
    axis.set(xticks=category, xticklabels=("Continuation", "Forward solves", "Adjoint work"), ylabel="Wall time (s)")
    axis.legend(fontsize=8)
    axis.grid(axis="y", alpha=.2)
    save(figure, directory, "optimization_runtime")


def fitted_report(fitted):
    if fitted is None:
        return ["The matched fitted-wall checks have not been supplied to this report."]
    lines = [
        "The linear and nonlinear optimization contours were reconstructed at density 0.5. "
        "Each contour was held fixed across refinements. Every fitted-wall result below uses "
        "Navier–Stokes at the common inlet velocity.",
        "",
        "| Geometry | Target spacing (mm) | Power (µW/m) | Resistance (Pa·s/m²) | Side flow (%) |",
        "|---|---:|---:|---:|---:|",
    ]
    for key, name in zip(("original", "stokes", "navier_stokes"), DESIGN_LABELS):
        for row in fitted["series"][key]:
            if row["state_passed"]:
                lines.append(f"| {name} | {row['mesh_target_mm']:g} | {1e6*row['mechanical_power_W_per_m']:.6f} "
                             f"| {row['hydraulic_resistance_pa_s_per_m2']:.4f} | {100*row['side_fraction']:.6f} |")
    comparison = fitted["comparison"]
    point = comparison.get("finest_paired")
    lines.append("")
    if point is None:
        lines.append("A matched finest pair was unavailable.")
        return lines
    power = point["stokes_minus_navier_stokes_power_W_per_m"] * 1e6
    side = point["stokes_minus_navier_stokes_side_flow_pp"]
    status = {
        "stokes_contour_lower_power": "The Stokes-optimized wall has lower pumping power under fitted Navier–Stokes CFD.",
        "navier_stokes_contour_lower_power": "The Navier–Stokes-optimized wall has lower pumping power under fitted Navier–Stokes CFD.",
        "unresolved": "The fitted power difference is unresolved under the paired numerical screen.",
        "incomplete": "The fitted power comparison is incomplete.",
    }
    lines.append(f"At {point['mesh_target_mm']:g} mm, the signed Stokes-minus–Navier–Stokes differences "
                 f"were {power:+.6f} µW/m in power and {side:+.6f} percentage points in side flow. "
                 + status[comparison["power_comparison_status"]])
    p_uncertainty = comparison.get("power_uncertainty_screen_W_per_m")
    s_uncertainty = comparison.get("side_flow_uncertainty_screen_pp")
    energy = comparison.get("combined_finest_raw_energy_residual_W_per_m")
    if p_uncertainty is not None and s_uncertainty is not None and energy is not None:
        lines.append(f"The paired uncertainty estimates were {p_uncertainty*1e6:.6f} µW/m and "
                     f"{s_uncertainty:.6f} percentage points. The sum of finest raw energy residual "
                     f"magnitudes was {energy*1e6:.6f} µW/m. Side-flow difference resolution under the "
                     f"screen was {'demonstrated' if comparison['side_flow_difference_resolved'] else 'unresolved'}.")
    lines.append("A power difference and a side-flow difference should be read together because the outlet "
                 "flow allocations can differ between the two optimized contours.")
    return lines


def write_report(output, result):
    cross, models = result["cross_evaluation"], result["optimization"]
    geometry = cross["geometry_comparison"]
    lines = [
        "# Does convective inertia change the optimized coolant-manifold junction?",
        "",
        "The Stokes comparison removes convective inertia from the momentum residual and its state Jacobian. "
        "The inlet velocity stays at the operating condition used for the Re = 50 Navier–Stokes model. "
        "The Brinkman term, zero-total-traction outlets, Taylor–Hood elements, anchored filter and editable cells "
        "stay common. The port objective includes kinetic-energy transport in both evaluators. The design problem "
        "remains nonlinear through the Brinkman density mapping and the performance functional.",
        "",
        "For each response, the adjoint solves the transpose of the free-DOF state Jacobian. The reduced gradient "
        "passes through the Brinkman interpolation, projection, anchored filter transpose and density transfer. "
        "Stokes has a linear state solve at each fixed density. Navier–Stokes requires nonlinear Newton iterations.",
        "",
        "## Cross-evaluation on the fixed-pipe Brinkman mesh",
        "",
        "Each physical model evaluates the original design and both optimized density fields. "
        "The common β = 8 side-flow floor is {:.6f}% and comes from the archived original Navier–Stokes state."
        .format(100 * cross["reference_side_fraction"]),
        "",
        "| Density design | Evaluator | Power (µW/m) | Resistance (Pa·s/m²) | Side flow (%) | Meets floor |",
        "|---|---|---:|---:|---:|---|",
    ]
    for row in cross["cases"]:
        design = DESIGN_LABELS[DESIGNS.index(row["design"])]
        model = MODEL_LABELS[MODELS.index(row["model"])]
        lines.append(f"| {design} | {model} | {1e6*row['power_w_per_m']:.6f} "
                     f"| {row['hydraulic_resistance_pa_s_per_m2']:.4f} | {100*row['side_fraction']:.6f} "
                     f"| {'Yes' if row['constraint_satisfied'] else 'No'} |")
    delta = cross["differences_by_evaluator"]["navier_stokes"]
    lines.extend((
        "",
        f"Evaluated under Navier–Stokes on the immersed mesh, the Stokes-optimized density changed power by "
        f"{delta['power_stokes_minus_navier_stokes_w_per_m']*1e6:+.6f} µW/m and side flow by "
        f"{delta['side_stokes_minus_navier_stokes_percentage_points']:+.6f} percentage points relative to the "
        "Navier–Stokes-optimized density.",
        "",
        "## Geometry",
        "",
        "![Original, Stokes-optimized and Navier–Stokes-optimized junctions](figures/junction_comparison.png)",
        "",
        "Across the 128 editable cells the two optima have a mean absolute raw density difference of "
        f"{geometry['mean_absolute_raw_density_difference']:.5f}, RMS difference of "
        f"{geometry['rms_raw_density_difference']:.5f}, and maximum difference of "
        f"{geometry['maximum_absolute_raw_density_difference']:.5f}. Their thresholded regions have a "
        f"symmetric-difference area of {geometry['symmetric_difference_area_mm2']:.6f} mm².",
        "",
        "| Geometry | Fluid area (mm²) | Maximum wall turn (degrees) |",
        "|---|---:|---:|",
    ))
    for key, label in zip(DESIGNS, DESIGN_LABELS):
        row = geometry["geometries"][key]
        lines.append(f"| {label} | {row['physical_area_mm2']:.6f} | {row['max_turn_degrees']:.6f} |")
    lines.extend(("", "## Fitted-wall Navier–Stokes comparison", ""))
    lines.extend(fitted_report(result["fitted_wall"]))
    lines.extend((
        "",
        "## Optimization cost and KKT",
        "",
        "![Continuation runtime, forward solves and adjoint work compared](figures/optimization_runtime.png)",
        "",
        "Continuation wall times include all reported stage work. Optimization work also depends on the number "
        "of function evaluations and any restoration. The archived and new run environments are retained in the JSON.",
        "",
        "| Model | Continuation wall time (s) | Forward (s) | Adjoint (s) | Unique MMA evaluations |",
        "|---|---:|---:|---:|---:|",
    ))
    for mode, name in zip(MODELS, MODEL_LABELS):
        cost = models[mode]["timing"]
        lines.append(f"| {name} | {cost['elapsed_seconds']:.3f} | {cost['forward_seconds']:.3f} "
                     f"| {cost['adjoint_seconds']:.3f} | {cost['mma_unique_evaluations']} |")
    linear = models["stokes"]["timing"]["elapsed_seconds"]
    nonlinear = models["navier_stokes"]["timing"]["elapsed_seconds"]
    if linear and nonlinear:
        lines.extend(("", f"The observed Stokes-to–Navier–Stokes continuation wall-time ratio was {linear/nonlinear:.3f}. "
                      "This includes differences in trajectory, derivative work, and iteration count."))
    lines.append("")
    for mode, name in zip(MODELS, MODEL_LABELS):
        if models[mode]["kkt"]:
            row = models[mode]["kkt"][-1]
            lines.append(f"- {name} endpoint audit {'passed' if row['verified'] else 'failed'}. "
                         f"Primal {row['primal_violation']:.3e}, dual {row['dual_violation']:.3e}, "
                         f"projected stationarity {row['stationarity_inf']:.3e}, and complementarity "
                         f"{row['complementarity_inf']:.3e}. Side-flow multiplier {row['multipliers'][0]:.5g}.")
    lines.extend(("", "## Figures", ""))
    for name in ("geometry_comparison", "junction_comparison", "optimized_walls_overlay", "stokes_density",
                 "navier_stokes_density", "density_difference", "velocity_comparison", "pressure_comparison",
                 "objective_history", "side_flow_history", "kkt_history", "optimization_runtime"):
        lines.append(f"- [{name.replace('_', ' ').capitalize()}]({FIGURE_DIR}/{name}.png)")
    content = "\n".join(lines) + "\n"
    (output / "nonlinear_comparison.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n",
                                                     encoding="utf-8")
    (output / "nonlinear_component.md").write_text(content, encoding="utf-8")
    reproduce = """
## Regenerate the comparisons

Run from the source tree with Python, NumPy, SciPy, NLopt, Shapely, Triangle and Matplotlib installed. Use one numerical-library thread for repeatable timing. The placeholders below refer to the archived pilot run, verification ZIP and the evidence directories.

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg
base=project4_fem.project43_topology.verification
python -B -u -m "$base.fixed_pipe_stokes_optimize" "$ns_pilot" stokes_run
python -B -u -m "$base.fixed_pipe_cross" "$ns_pilot" stokes_run/stage_03_beta_8/selected/fields.npz cross_evaluation
python -B -m "$base.plot_nonlinear_comparison" cross_evaluation stokes_run "$ns_pilot" "$ns_kkt_archive" report --fitted "$fitted_report"
```

See the fitted-wall run logs and source for the three Sol mesh evaluations. The six immersed field files contain actual P1 pressure and P2 velocity degrees of freedom. The velocity plots render the nodal speed on midpoint-subdivided triangles. The field files retain the vector components for exact P2 interpolation.
"""
    (output / "README.md").write_text(content + reproduce, encoding="utf-8")


def build(root, stokes, ns, ns_kkt, fitted, output):
    output.mkdir(parents=True, exist_ok=True)
    figure_dir = output / FIGURE_DIR
    figure_dir.mkdir(exist_ok=True)
    cross = read_json(root / "cross_evaluation.json")
    geometries = [shape(read_json(root / f"{design}.geometry.json")["geometry"]) for design in DESIGNS]
    fields = {}
    for case in cross["cases"]:
        with np.load(root / case["fields_file"], allow_pickle=False) as data:
            fields[(case["design"], case["model"])] = {key: data[key] for key in data.files}
    models = {"stokes": model_data(stokes, None), "navier_stokes": model_data(ns, ns_kkt)}
    geometry_figures(geometries, figure_dir)
    density_figures(fields, geometries, figure_dir)
    field_figures(fields, geometries, figure_dir)
    history_figures(models, figure_dir)
    result = {"cross_evaluation": cross, "optimization": models,
              "fitted_wall": read_json(fitted) if fitted is not None else None}
    write_report(output, result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cross", type=Path)
    parser.add_argument("stokes", type=Path)
    parser.add_argument("navier_stokes", type=Path)
    parser.add_argument("ns_kkt", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--fitted", type=Path)
    args = parser.parse_args()
    build(args.cross, args.stokes, args.navier_stokes, args.ns_kkt, args.fitted, args.output)


if __name__ == "__main__":
    main()
