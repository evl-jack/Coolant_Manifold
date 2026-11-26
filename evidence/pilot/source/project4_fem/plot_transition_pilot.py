"""Presentation figures based on the measured pilot and marked analytical proxies."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.tri import Triangulation
import numpy as np

from .fillet import BASELINE, FilletMap
from .plot_recovery_audit import BLUE, BROWN, INK, TEAL, geometry_figure
from .recovery_audit import GEOMETRY, theoretical_authority
from .transition_pilot import LOWER, UPPER, TransitionMap


COLORS = (BLUE, TEAL)
LABELS = ("Stokes", "Re 50")


def save(fig, path):
    fig.savefig(path.with_suffix(".png"), dpi=210, facecolor="white")
    fig.savefig(path.with_suffix(".pdf"), facecolor="white")
    plt.close(fig)


def axes(ax):
    ax.axhspan(100*LOWER, 100*UPPER, color="#dcece2")
    ax.text(.02, 100*(LOWER+UPPER)/2, "19.5%–20.5% target", transform=ax.get_yaxis_transform(),
            va="center", color=INK, fontsize=9)
    ax.grid(axis="y", linewidth=.6, color="#e6ebef")
    ax.set_ylabel("Side share of inlet flow (%)")
    ax.spines[["top", "right"]].set_visible(False)


def point(ax, x, row, color, shift=0.):
    if not row["converged"]:
        return
    y = 100*row["side_share"]
    marker = "D" if row["case"] == "baseline" else "o"
    ax.scatter(x+shift, y, marker=marker, s=58, facecolors=color if row["area_feasible"] else "white",
               edgecolors=color, linewidth=1.4, zorder=3)
    if not row["pressure_feasible"]:
        ax.scatter(x+shift, y, s=105, facecolors="none", edgecolors=BROWN, linewidth=1.1, zorder=4)


def legend(ax):
    handles = [Line2D([], [], marker="o", linestyle="none", color=c, label=l) for c, l in zip(COLORS, LABELS)]
    handles += [Line2D([], [], marker="D", linestyle="none", color=INK, label="Baseline"),
                Line2D([], [], marker="o", markerfacecolor="white", linestyle="none", color=INK, label="Area violation"),
                Line2D([], [], marker="o", markerfacecolor="white", linestyle="none", color=BROWN, label="Pressure violation ring")]
    ax.legend(handles=handles, ncol=3, loc="lower center", bbox_to_anchor=(.5, .19),
              frameon=False, fontsize=9)


def measured(models, output):
    fig, ax = plt.subplots(figsize=(11.6, 5), layout="constrained")
    axes(ax)
    short = ["Baseline", "Narrow 100%", "Narrow 75%", "Narrow 50%", "Narrow 25%", "Widen 50%", "Widen 100%",
             "Curvature", "Archived probe", "Base + 50%", "Base + 100%", "Probe + 50%", "Probe + 100%"]
    for color, shift, cases in zip(COLORS, (-.12, .12), models):
        for i, row in enumerate(cases):
            point(ax, i, row, color, shift)
    ax.set(xticks=range(len(short)), xticklabels=short, title="Measured 26-control FEM pilot",
           xlabel="Combined cases use area compensation", ylim=(18.8, 29.2))
    ax.tick_params(axis="x", labelrotation=35)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    legend(ax)
    save(fig, output/"measured-side-flow")


def prediction(models, output):
    fig, panels = plt.subplots(1, 2, figsize=(10.8, 4.55), sharey=True, layout="constrained")
    for ax, name, color, cases in zip(panels, LABELS, COLORS, models):
        axes(ax)
        sweeps = sorted((r for r in cases if r["mode"] == "width sweep" and r["converged"]), key=lambda row: row["transition_dial"])
        baseline = next(r for r in sweeps if r["transition_dial"] == 0.)
        proxy = theoretical_authority(FilletMap(GEOMETRY), BASELINE, baseline["side_share"])
        peak = proxy["max_half_width_excursion_mm"]
        ax.plot([peak*r["dial"] for r in proxy["cases"]], [100*r["side_share"] for r in proxy["cases"]],
                "--", color=BROWN, label="Friction-only estimate")
        ax.plot([peak*r["transition_dial"] for r in sweeps], [100*r["side_share"] for r in sweeps],
                "o-", color=color, label="Measured FEM")
        ax.set(title=name, xlabel="Signed peak half-width displacement (mm)", ylim=(18.9, 29.0))
        ax.legend(frameon=False, loc="lower right", bbox_to_anchor=(1., .19), fontsize=9)
    save(fig, output/"prediction-versus-fem")


def area(models, output):
    fig, ax = plt.subplots(figsize=(8.8, 5.1), layout="constrained")
    axes(ax)
    ax.axvline(0., color="#748797", linestyle="--", linewidth=1.)
    for color, cases in zip(COLORS, models):
        for row in cases:
            point(ax, row["physical_fluid_area_change_mm2"], row, color)
    ax.set(title="Area and side-flow response", xlabel="Physical fluid-area change from matched baseline (mm²)", ylim=(18.8, 29.2))
    legend(ax)
    save(fig, output/"area-side-flow")


def velocity(models, output):
    with np.load(output/"velocity-re-50.npz") as data:
        case = str(data["case"])
        selected = next(r for r in models[1] if r["case"] == case)
        points = 2*data["velocity_points"]+[-1., 3.]
        a, b, c, ab, bc, ca = data["velocity_dofs"].T
        elements = np.vstack([np.c_[a, ab, ca], np.c_[ab, b, bc], np.c_[ca, bc, c], np.c_[ab, bc, ca]])
        n = len(points)
        speed = np.hypot(data["state"][:n], data["state"][n:2*n])
    mapping = TransitionMap(GEOMETRY)
    y = np.linspace(5., 8., 600)
    width = mapping.basis(y)[0]+mapping.offset(y, np.array(selected["design"]))
    fig, ax = plt.subplots(figsize=(10.5, 4.65), layout="constrained")
    field = ax.tripcolor(Triangulation(*points.T, elements), speed, shading="gouraud", cmap="viridis", rasterized=True)
    for sign in (-1, 1):
        ax.plot(10+sign*width, y, color="white", lw=1, ls="--", label="Physical wall contour" if sign == -1 else None)
    fig.colorbar(field, ax=ax, label="Speed / mean inlet speed", shrink=.9)
    ax.set(aspect="equal", xlabel="x (mm)", ylabel="y (mm)", title="Converged Re 50 pilot velocity field")
    ax.text(.01,.96,case,transform=ax.transAxes,va="top",color=INK)
    ax.legend(loc="upper right", facecolor=INK, labelcolor="white", fontsize=9)
    save(fig, output/"pilot-velocity-re-50")


def plot(output):
    plt.rcParams.update({"font.size":11, "axes.titlesize":14, "figure.facecolor":"white", "text.color":INK})
    models = [json.loads((output/f"cases-re-{re}.json").read_text())["rows"] for re in (0, 50)]
    geometry_figure(output, pilot=True)
    measured(models, output)
    prediction(models, output)
    area(models, output)
    velocity(models, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    plot(parser.parse_args().output)
