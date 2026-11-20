"""Plot the saved shape candidate on the explicit-wall mesh."""

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import numpy as np

from .shape_search import CASE, NAMES, flow_metrics, manifold_grid_axes, solve_walls


def main():
    folder = Path(__file__).resolve().parent / "results"
    with (folder / "shape_best.csv").open(newline="") as file:
        best = next(csv.DictReader(file))
    controls = [float(best[name]) for name in NAMES]
    spacing = .00025
    mask, result = solve_walls(controls, spacing)
    target = CASE["side_to_main_ratio"]
    metrics = flow_metrics(result, target)
    x, y = (edges * 1e3 for edges in manifold_grid_axes(spacing))

    styles = {"wall": ("#263643", "Wall"), "inlet": ("#73b3de", "Inlet"),
              "main_outlet": ("#e39346", "Main outlet"), "side_outlet": ("#bf80d5", "Side outlet")}
    segments = {tag: [] for tag in styles}
    for face in result["boundary"]:
        nx, ny = face["normal"]
        half = face["length"] / 2 * np.array([-ny, nx])
        segments[face["tag"]].append(1e3 * np.array([face["center"] - half, face["center"] + half]))

    fig, axes = plt.subplots(2, 1, figsize=(9, 8), layout="constrained")
    for ax, values, title, units, colors in (
        (axes[0], np.hypot(result["u"], result["v"]), "Velocity magnitude and direction", "Speed (m/s)", "viridis"),
        (axes[1], result["p"], "Static gauge pressure", "Pressure (Pa)", "coolwarm"),
    ):
        field = np.full(mask.shape, np.nan)
        field[mask] = values
        colors_plot = ax.pcolormesh(x, y, field.T, cmap=colors, shading="flat")
        fig.colorbar(colors_plot, ax=ax, label=units, pad=.025, shrink=.85)
        for tag, (color, _) in styles.items():
            ax.add_collection(LineCollection(segments[tag], colors=[color],
                                            linewidths=1.0 if tag == "wall" else 2.5, zorder=3))
        ax.set(xlabel="x (mm)", ylabel="y (mm)", xlim=(-.3, x[-1] + .3),
               ylim=(1.5, y[-1] + .5), aspect="equal", title=title)

    sample = np.zeros(mask.shape, dtype=bool)
    sample[1::3, 1::3] = True
    sample = sample[mask]
    xx, yy = result["centers"][sample].T * 1e3
    arrows = axes[0].quiver(xx, yy, result["u"][sample], result["v"][sample],
                            angles="xy", scale_units="xy", scale=7, color="white",
                            edgecolors="#20303a", linewidths=.3, zorder=4)
    axes[0].quiverkey(arrows, 16, 8.5, 2, "2 m/s", coordinates="data", color="#263643")
    axes[1].legend(handles=[Line2D([], [], color=color, linewidth=2, label=label)
                           for color, label in styles.values()],
                   loc="upper left", ncol=2, frameon=False)
    status = "converged" if result["converged"] else "unconverged"
    fig.suptitle(f"Manifold shape candidate · explicit no-slip walls · {spacing * 1e3:g} mm mesh\n"
                 f"Side/main flow ratio {metrics['side_to_main_ratio']:.3f} of {target:.3f} target · "
                 f"{status} in {result['iterations']} iterations")
    image_path = folder / "shape_flow.png"
    fig.savefig(image_path, dpi=180)
    plt.close(fig)
    print(f"Saved {image_path}")
    print(f"Side share {metrics['side_split']:.2%}, side/main {metrics['side_to_main_ratio']:.4f}, {status}")


if __name__ == "__main__":
    main()
