"""Solve the saved optimized density and create its flow/pressure plot."""
import argparse
import csv
from pathlib import Path

import numpy as np

try:
    from .solve import solve
    from .case import CASE
    from .optimize import flow_metrics, outlet_violation
except ImportError:
    from solve import solve
    from case import CASE
    from optimize import flow_metrics, outlet_violation

from mesh import build_mesh, manifold_grid_axes


def save_plot(result, density, path, ratio, target):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D

    centers = result["centers"]
    x, y = np.unique(centers[:, 0]), np.unique(centers[:, 1])
    shape = (len(x), len(y))
    fields = {name: result[name].reshape(shape).T for name in ("u", "v", "p")}
    gamma = np.asarray(density).reshape(shape).T
    speed = np.hypot(fields["u"], fields["v"])
    xm, ym = x * 1e3, y * 1e3
    x_edges, y_edges = (axis * 1e3 for axis in manifold_grid_axes(CASE["spacing_m"]))
    left, right = CASE["main_channel"]["x"]
    fixed_run = CASE["fixed_main_run_m"]
    fixed_edges = np.array([left + fixed_run, right - fixed_run]) * 1e3
    fixed_side_start = (CASE["side_branch"]["y"][1] - CASE["fixed_side_run_m"]) * 1e3
    reference = []
    for face in build_mesh(CASE["spacing_m"])["boundary_faces"]:
        if face["tag"] == "wall":
            nx, ny = face["normal"]
            half = face["length"] / 2 * np.array([-ny, nx])
            reference.append(np.array([face["center"] - half, face["center"] + half]) * 1e3)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), constrained_layout=True)
    for ax, values, title, cmap in (
        (axes[0], speed, "Speed (m/s)", "viridis"),
        (axes[1], fields["p"], "Gauge pressure (Pa)", "coolwarm"),
    ):
        image = ax.pcolormesh(x_edges, y_edges, values, shading="flat", cmap=cmap)
        fig.colorbar(image, ax=ax, label=title)
        ax.add_collection(LineCollection(reference, colors="white", linestyles="dashed", linewidths=1.1))
        ax.contour(xm, ym, gamma, levels=[.5], colors="black", linewidths=1.2)
        ax.vlines(fixed_edges, y_edges[0], y_edges[-1], colors=".65", linestyles="dotted", linewidth=1)
        ax.hlines(fixed_side_start, x_edges[0], x_edges[-1], colors=".65", linestyles="dotted", linewidth=1)
        ax.set(title=title, xlabel="x (mm)", ylabel="y (mm)", aspect="equal")
    axes[1].legend(handles=[Line2D([], [], color="white", linestyle="--", label="Starting T walls"),
                           Line2D([], [], color="black", label="Final density contour"),
                           Line2D([], [], color=".4", linestyle=":", label="Fixed section limits")],
                   loc="upper right", facecolor=".65", framealpha=1)

    # Sparse vectors show direction on this nonuniform grid.
    xx, yy = np.meshgrid(xm, ym, indexing="ij")
    take = np.zeros(shape, dtype=bool)
    take[::2, ::2] = True
    take &= gamma.T >= .5
    arrows = axes[0].quiver(xx[take], yy[take], fields["u"].T[take] * 1e3,
                            fields["v"].T[take] * 1e3, color="white",
                            angles="xy", scale_units="xy", scale=2500)
    axes[0].quiverkey(arrows, .88, .94, 2000, "2 m/s", coordinates="axes",
                      color="black", labelcolor="black")
    solver_status = "converged" if result["converged"] else "unconverged"
    fig.suptitle(f"2D RANS flow ({solver_status} after {result['iterations']} iterations)\n"
                 f"Side / main outlet flow {ratio:.3f} · target {target:.3f}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main(density_path=None, max_iterations=1000, show=False):
    folder = Path(__file__).resolve().parent / "results"
    target = CASE["side_to_main_ratio"]
    tolerance = CASE["outlet_balance_tolerance"]
    model = "B"
    history_path = folder / "optimization.csv"
    if density_path is None and history_path.exists():
        with history_path.open(newline="") as file:
            rows = list(csv.DictReader(file))
        if rows:
            target = float(rows[-1].get("target_side_to_main_ratio", target))
            tolerance = float(rows[-1].get("outlet_balance_tolerance", tolerance))
            model = rows[-1].get("model", model)
    density_path = Path(density_path) if density_path else folder / "optimized_density.npy"
    density = np.load(density_path)
    result = solve(max_iterations=max_iterations, model=model, design_density=density,
                   alpha_max=1e6, convergence_tolerance=5e-4)
    metrics = flow_metrics(result, target)
    ratio = metrics["side_to_main_ratio"]
    feasible = result["converged"] and outlet_violation(metrics, tolerance) == 0

    stem = folder / "optimized_flow"
    np.savez_compressed(stem.with_suffix(".npz"), density=density,
                        converged=result["converged"], iterations=result["iterations"],
                        target_side_to_main_ratio=target, outlet_balance_tolerance=tolerance,
                        feasible=feasible, **metrics,
                        **{name: result[name] for name in ("centers", "areas", "u", "v", "p", "k", "omega", "nut")})
    save_plot(result, density, stem.with_suffix(".png"), ratio, target)
    print(f"Flow fields saved to {stem.with_suffix('.npz')}")
    print(f"Plot saved to {stem.with_suffix('.png')}")
    print(f"Converged: {result['converged']} after {result['iterations']} iterations")
    print(f"Side / main flow: {ratio:.4f} (target {target:.4f}), side share: {metrics['side_split']:.2%}")
    print(f"Outlet flow (kg/s per m depth): side {metrics['side_mass_flow']:.4f}, main {metrics['main_mass_flow']:.4f}")
    if show:
        import os
        os.startfile(stem.with_suffix(".png"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--density", help="density .npy file; defaults to results/optimized_density.npy")
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--show", action="store_true", help="open the saved plot")
    args = parser.parse_args()
    main(args.density, args.iterations, args.show)
