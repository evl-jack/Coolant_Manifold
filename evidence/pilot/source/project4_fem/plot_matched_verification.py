"""Plot measured FEM states from the matched wall-treatment study."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def plot(directory):
    data = json.loads((directory/"matched-verification.json").read_text())["rows"]
    colors = {"Baseline": "#315c97", "Candidate": "#087d75"}
    styles = {"Brinkman": ("-", "o"), "Fitted": ("--", "s")}
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.8), sharey=True, layout="constrained")
    fig.suptitle("Matched FEM refinement and wall treatment", fontsize=15)
    for ax, reynolds, label in zip(axes, (0, 50), ("Stokes", "Re 50")):
        for method, (line, marker) in styles.items():
            for shape, color in colors.items():
                rows = [row for row in data if row["reynolds"] == reynolds
                        and row["method"] == method and row["shape"] == shape]
                rows.sort(key=lambda row: row["level"])
                ax.plot([r["level"] for r in rows],
                        [100*r["side_share"] if r["converged"] else float("nan") for r in rows],
                        color=color, linestyle=line, marker=marker, lw=1.7, ms=5)
        ax.axhspan(19.5, 20.5, color="#dcece2")
        ax.text(.03, 20., "19.5%–20.5% target", transform=ax.get_yaxis_transform(),
                va="center", fontsize=9, color="#244c36")
        ax.grid(axis="y", color="#e8ecf0", lw=.7)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set(title=label, xticks=[1, 2, 3], xticklabels=["2,500", "10,000", "22,500"],
               xlabel="Mesh triangles", xlim=(.82, 3.18), ylim=(19, 29))
    axes[0].set_ylabel("Side share of inlet flow (%)")
    keys = [Line2D([], [], color=color, lw=2, label=shape) for shape, color in colors.items()]
    keys += [Line2D([], [], color="#31404b", linestyle=line, marker=marker, label=method)
             for method, (line, marker) in styles.items()]
    fig.legend(handles=keys, loc="outside lower center", ncol=4, frameon=False)
    fig.savefig(directory/"matched-side-flow.png", dpi=240, facecolor="white")
    fig.savefig(directory/"matched-side-flow.pdf", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    plot(parser.parse_args().directory)
