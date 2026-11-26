"""Render archived and proposed geometry diagnostics for the decision report."""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np

from .fillet import FilletMap
from .recovery_audit import GEOMETRY, TRANSITION_LOCATIONS, transition


INK, BLUE, TEAL, BROWN = "#162b3d", "#305b99", "#087f78", "#99482f"


def geometry_figure(output, pilot=False):
    m = FilletMap(GEOMETRY)
    fig, ax = plt.subplots(figsize=(9.2, 5.4), layout="constrained")
    y = np.linspace(5., 6.5, 500)
    _, clearance, amplitude = m.parameters(y)
    baseline = .5+clearance
    yt = np.linspace(6.5, 8., 400)
    throat, added = transition(yt)
    knots = m.spline.t
    arc = np.array([np.mean(knots[i+1:i+4]) for i in range(20)])
    old_nodes = 5+1.5*(1-np.cos(.5*np.pi*arc))
    centers = .5+m.parameters(old_nodes)[1]
    new_nodes = transition(TRANSITION_LOCATIONS)[0]
    for sign in (-1, 1):
        left = sign < 0
        ax.fill_betweenx(y, 10+sign*(baseline-1.14*amplitude), 10+sign*(baseline+1.14*amplitude),
                         color="#a7c4e7", alpha=.8, label="Existing fillet density strips" if left else None)
        ax.fill_betweenx(yt, 10+sign*(throat-added), 10+sign*(throat+added),
                         color="#9bdbd0", alpha=.85,
                         label=("Added physical transition envelope" if pilot else "Proposed physical transition envelope") if left else None)
        ax.plot(10+sign*baseline, y, color=INK, lw=1.3)
        ax.plot(10+sign*throat, yt, color=INK, lw=1.3)
        ax.plot([10+sign*.625]*2, [8., 8.45], color=INK, lw=1.3)
        ax.hlines(5., *(sorted((10+sign*2, 10+sign*2.5))), color=INK, lw=1.3)
        ax.scatter(10+sign*centers, old_nodes, s=13, color=BLUE, zorder=3, label="20 mirrored fillet controls" if left else None)
        ax.scatter(10+sign*new_nodes, TRANSITION_LOCATIONS, s=31, marker="s", facecolor="white", edgecolor=TEAL,
                   zorder=4, label=("6 added mirrored controls" if pilot else "6 proposed mirrored controls") if left else None)
    ax.axhline(6.5, color="#a7b4bf", lw=.8, ls="--")
    ax.axhline(8., color="#a7b4bf", lw=.8, ls="--")
    ax.text(7.68, 8.05, "Fixed outlet guide", color=INK, fontsize=9)
    ax.set(xlabel="x (mm)", ylabel="y (mm)", xlim=(7.45, 12.55), ylim=(4.85, 8.45), aspect="equal",
           title=("26-control pilot envelope and fixed joins" if pilot else "Proposed transition controls and existing fillet envelope"))
    ax.legend(loc="lower center", bbox_to_anchor=(.5, -.26), ncol=2, frameon=False, fontsize=9)
    fig.savefig(output/"design-envelopes.png", dpi=175)
    if pilot:
        fig.savefig(output/"design-envelopes.pdf")
    plt.close(fig)


def authority_figure(data, output):
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.6, 4.25), sharey=True, layout="constrained")
    for ax in (a, b):
        ax.axhspan(19.5, 20.5, color="#dcece2", label="Flow target band")
        ax.grid(axis="y", alpha=.2)
        ax.set_ylim(18.9, 29.)
    for name, label, c in (("re_0", "Stokes", BLUE), ("re_50", "Re 50", TEAL)):
        points = data[name]["amplitude"]
        a.plot([r["amplitude_mm"] for r in points], [100*r["side_share"] for r in points], "o-", color=c, label=label)
    a.set(title="Existing fillet · measured FEM response", xlabel="Signed peak half-width change (mm)",
          ylabel="Side share of inlet flow (%)")
    a.legend(frameon=False, loc="lower right")
    proxy = data["transition_proxy"]
    peak = proxy["max_half_width_excursion_mm"]
    b.plot([peak*r["dial"] for r in proxy["cases"]], [100*r["side_share"] for r in proxy["cases"]],
           "o--", markerfacecolor="white", color=BROWN, label="Friction-only estimate")
    b.plot(0., 100*proxy["start_share"], marker="D", color=INK, label="Starting explicit probe")
    b.set(title="Proposed transition · analytical estimate", xlabel="Signed peak half-width change (mm)")
    b.legend(frameon=False, loc="lower right")
    fig.savefig(output/"side-flow-authority.png", dpi=175)
    plt.close(fig)


def velocity_figure(state_file, output):
    with np.load(state_file) as archive:
        points = 2*archive["velocity_points"] + [-1., 3.]
        a, b, c, ab, bc, ca = archive["velocity_dofs"].T
        elements = np.vstack([np.c_[a, ab, ca], np.c_[ab, b, bc], np.c_[ca, bc, c], np.c_[ab, bc, ca]])
        n = len(points)
        speed = np.hypot(archive["state"][:n], archive["state"][n:2*n])
    fig, ax = plt.subplots(figsize=(10.5, 4.7), layout="constrained")
    field = ax.tripcolor(Triangulation(*points.T, elements), speed, shading="gouraud", cmap="viridis")
    fig.colorbar(field, ax=ax, label="Speed / mean inlet speed", shrink=.9)
    ax.set(aspect="equal", xlabel="x (mm)", ylabel="y (mm)", title="Frozen FEM reference · Re 50 velocity magnitude")
    fig.savefig(output/"baseline-velocity.png", dpi=175)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    args = parser.parse_args()
    geometry_figure(args.output)
    authority_figure(json.loads((args.output/"audit.json").read_text()), args.output)
    velocity_figure(args.baseline, args.output)
