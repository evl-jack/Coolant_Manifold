"""Differentiable, leakage-independent connectivity certificate for the T-junction.

For each outlet, freeze inlet potential at one and that outlet at zero. The
other outlet floats. Minimizing the graph energy gives the effective port
conductance. Each edge has a positive numerical floor and a smooth term that
is exactly zero unless both of its ends exceed the conservative fluid cutoff.
If the positive-gated graph is disconnected, the inlet component has a cut
whose current is bounded above by the total floor budget. A conductance above
that budget proves that a face-connected fluid path exists. This graph does
not use the Brinkman velocity, pressure, resistance, or outlet flux.
"""

import argparse
from dataclasses import dataclass, replace
import json
import os
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve


@dataclass(frozen=True)
class ConnectivityResult:
    conductance: np.ndarray
    projected_jacobian: np.ndarray
    constraints: np.ndarray
    constraint_projected_jacobian: np.ndarray
    hard_paths: tuple[bool, bool]
    height_jacobian: np.ndarray
    constraint_height_jacobian: np.ndarray


def _gate(density, cutoff, width):
    t = np.clip((density - cutoff) / width, 0.0, 1.0)
    smooth = t**3 * (10.0 + t * (-15.0 + 6.0 * t))
    derivative = 30.0 * t**2 * (1.0 - t)**2 / width
    return smooth, derivative


@dataclass(frozen=True)
class RectangularBaffleClearance:
    """Smooth edge gates for a sharp rectangular wall-attached baffle.

    A whole edge has positive gate only if it clears a lateral side, the
    base, or the tip by more than the configured margin. The lateral/base
    clearances are constant. The tip clearance decreases by one mm per mm
    of height. Their smooth union never opens an edge through the rectangle.
    """

    x_min_mm: float = 13.0
    x_max_mm: float = 14.0
    wall_y_mm: float = 0.0
    direction: int = 1
    margin_mm: float = 0.025
    transition_mm: float = 0.125
    max_height_mm: float = 3.75

    def __post_init__(self):
        values = [self.x_min_mm, self.x_max_mm, self.wall_y_mm, self.margin_mm,
                  self.transition_mm, self.max_height_mm]
        if (not np.all(np.isfinite(values)) or self.x_min_mm >= self.x_max_mm
                or self.direction not in (-1, 1) or min(self.margin_mm, self.transition_mm, self.max_height_mm) <= 0):
            raise ValueError("Invalid baffle bounds, orientation or clearance window")

    def fixed_clearances(self, segments):
        x = segments[..., 0]
        s = self.direction * (segments[..., 1] - self.wall_y_mm)
        lateral = np.maximum(self.x_min_mm - np.max(x, axis=1), np.min(x, axis=1) - self.x_max_mm)
        lateral_gate, _ = _gate(lateral, self.margin_mm, self.transition_mm)
        base_gate, _ = _gate(-np.max(s, axis=1), self.margin_mm, self.transition_mm)
        return (1.0 - lateral_gate) * (1.0 - base_gate), np.min(s, axis=1)

    def evaluate(self, edge_segments_mm, height_mm):
        return self.from_fixed_clearances(self.fixed_clearances(edge_segments_mm), height_mm)

    def from_fixed_clearances(self, constant, height_mm):
        if not np.isfinite(height_mm) or not -1e-12 <= height_mm <= self.max_height_mm + 1e-12:
            raise ValueError(f"Baffle height must be in [0, {self.max_height_mm:g}] mm")
        closed_except_tip, tip_coordinate = constant
        tip_gate, tip_derivative = _gate(tip_coordinate - height_mm, self.margin_mm, self.transition_mm)
        gate = (1.0 - closed_except_tip) + closed_except_tip * tip_gate
        return gate, -closed_except_tip * tip_derivative


class ConnectivityGraph:
    """Face-adjacent sample graph with distinct fixed-fluid inlet and outlet nodes.

    Seams join a port to its design boundary only where the center's normal
    projection lies in the interior of the exact fixed-fluid attachment.
    An expanded envelope must use the same fixed-strip and constant-normal
    extension conventions in its reconstruction for this proof to apply.
    """

    def __init__(self, grid, *, threshold=0.5, margin=0.025, transition_width=0.15,
                 floor_budget=1e-4, minimum_conductance=1e-3, baffle=None):
        if (not np.all(np.isfinite([threshold, margin, transition_width, floor_budget,
                                   minimum_conductance]))
                or not 0 < threshold < threshold + margin < threshold + margin + transition_width <= 1
                or not 0 < floor_budget < minimum_conductance):
            raise ValueError("Invalid gate interval or certificate threshold")
        bounds = (grid.x_edges_mm[0], grid.x_edges_mm[-1], grid.y_edges_mm[0], grid.y_edges_mm[-1])
        west, east, bottom, top = bounds
        if not (0 < west < 9.5 < 10.5 < east < 20 and bottom <= 0 and 4 <= top < 7):
            raise ValueError("The rectangular design envelope must retain all three fixed-fluid seams")
        self.size = grid.size
        self.ports = tuple(range(self.size, self.size + 3))
        self.bounds_mm = bounds
        self.threshold = float(threshold)
        self.cutoff = float(threshold + margin)
        self.transition_width = float(transition_width)
        self.floor_budget = float(floor_budget)
        self.minimum_conductance = float(minimum_conductance)
        self.baffle = baffle
        if baffle is not None:
            far_y = baffle.wall_y_mm + baffle.direction * baffle.max_height_mm
            if not (west < baffle.x_min_mm < baffle.x_max_mm < east
                    and bottom <= baffle.wall_y_mm <= top and bottom <= far_y <= top):
                raise ValueError("The clearance proof requires the whole baffle within the design envelope")
        cells = np.arange(self.size).reshape(grid.shape)
        x, y = grid.centers_mm.T
        pairs = [np.column_stack((cells[:, :-1].ravel(), cells[:, 1:].ravel())),
                 np.column_stack((cells[:-1, :].ravel(), cells[1:, :].ravel()))]
        seams = (cells[:, 0][(y[cells[:, 0]] > 0) & (y[cells[:, 0]] < 4)],
                 cells[:, -1][(y[cells[:, -1]] > 0) & (y[cells[:, -1]] < 4)],
                 cells[-1, :][(x[cells[-1, :]] > 9.5) & (x[cells[-1, :]] < 10.5)])
        counts = [len(seam) for seam in seams]
        if min(counts) == 0:
            raise ValueError("At least one exact fixed-fluid center-to-face attachment is required per port")
        edge_weights = [np.ones(len(edges)) for edges in pairs]
        for port, seam in zip(self.ports, seams):
            pairs.append(np.column_stack((seam, np.full(len(seam), port))))
            edge_weights.append(np.full(len(seam), 2.0))
        self.seam_counts = tuple(counts)
        self.ends = np.vstack(pairs)
        self.weights = np.concatenate(edge_weights)
        if baffle is not None:
            segments = [grid.centers_mm[edges] for edges in pairs[:2]]
            for axis, location, seam in zip((0, 0, 1), (west, east, top), seams):
                centers = grid.centers_mm[seam]
                exterior = centers.copy()
                exterior[:, axis] = location
                segments.append(np.stack((centers, exterior), axis=1))
            self._baffle_clearance = baffle.fixed_clearances(np.vstack(segments))
        self.edge_floor = self.floor_budget / len(self.ends)
        rows = np.repeat(np.arange(len(self.ends)), 2)
        values = np.tile([-1., 1.], len(self.ends))
        self.incidence = coo_matrix((values, (rows, self.ends.ravel())),
                                    shape=(len(self.ends), self.size + 3)).tocsr()

    def _solve(self, gates, derivative, baffle_height_mm=0.0):
        left, right = self.ends.T
        if self.baffle is None:
            if baffle_height_mm != 0.0:
                raise ValueError("Baffle height provided for a graph without a baffle")
            edge_gate = np.ones(len(self.ends))
            edge_gate_derivative = np.zeros(len(self.ends))
        else:
            edge_gate, edge_gate_derivative = self.baffle.from_fixed_clearances(self._baffle_clearance, baffle_height_mm)
        active_edges = (gates[left] > 0.) & (gates[right] > 0.) & (edge_gate > 0.)
        adjacency = coo_matrix((np.ones(np.count_nonzero(active_edges)),
                                (left[active_edges], right[active_edges])),
                               shape=(self.size + 3, self.size + 3)).tocsr()
        _, components = connected_components(adjacency, directed=False)
        source, *sinks = self.ports
        hard_paths = tuple(bool(components[source] == components[sink]) for sink in sinks)
        gated_weights = self.weights * edge_gate
        gate_product = gates[left] * gates[right]
        edge_g = self.edge_floor + gated_weights * gate_product
        edge_g_height_derivative = self.weights * gate_product * edge_gate_derivative
        laplacian = (self.incidence.T @ self.incidence.multiply(edge_g[:, None])).tocsc()
        scores = np.empty(2)
        jac = np.zeros((2, self.size + 3))
        height_jac = np.empty(2)
        for i, sink in enumerate(sinks):
            free = np.setdiff1d(np.arange(self.size + 3), [source, sink])
            phi = np.zeros(self.size + 3)
            phi[source] = 1.0
            phi[free] = spsolve(laplacian[free][:, free], -laplacian[free, source].toarray().ravel())
            if not np.all(np.isfinite(phi)) or np.any(phi < -1e-8) or np.any(phi > 1.0 + 1e-8):
                raise RuntimeError("Graph Dirichlet solve violated the discrete maximum principle")
            drop_squared = (self.incidence @ phi)**2
            scores[i] = edge_g @ drop_squared
            height_jac[i] = edge_g_height_derivative @ drop_squared
            np.add.at(jac[i], left, drop_squared * gated_weights * derivative[left] * gates[right])
            np.add.at(jac[i], right, drop_squared * gated_weights * gates[left] * derivative[right])
        return scores, jac[:, :self.size], hard_paths, height_jac

    def evaluate(self, projected, baffle_height_mm=0.0):
        density = np.asarray(projected, dtype=float)
        if (density.shape != (self.size,) or not np.all(np.isfinite(density))
                or np.any(density < -1e-12) or np.any(density > 1.0 + 1e-12)):
            raise ValueError("Expected one finite projected density in [0, 1] per design cell")
        gates, derivative = _gate(density, self.cutoff, self.transition_width)
        scores, jac, hard, height_jac = self._solve(np.r_[gates, 1., 1., 1.], np.r_[derivative, 0., 0., 0.],
                                                   baffle_height_mm)
        constraints = self.minimum_conductance - scores
        if any(g <= 0 and not connected for g, connected in zip(constraints, hard)):
            raise RuntimeError("A computed conductance pass contradicts the graph certificate")
        return ConnectivityResult(scores, jac, constraints, -jac, hard, height_jac, -height_jac)

    def evaluate_density(self, pipeline, state, baffle_height_mm=0.0):
        """Return the certificate and its raw-density constraint Jacobian."""
        if pipeline.threshold != self.threshold:
            raise ValueError("Density/reconstruction and connectivity thresholds disagree")
        result = self.evaluate(state.projected, baffle_height_mm)
        raw_jacobian = np.asarray([pipeline.filter_pullback(row * state.projection_derivative)
                                   for row in result.constraint_projected_jacobian])
        return result, raw_jacobian


class ConnectivityDiagnosticEvaluator:
    """Opt-in full-density adapter for the existing diagnostic optimizer."""

    def __init__(self, hydraulics, graph=None):
        if hydraulics.reduced is not None:
            raise ValueError("Connectivity diagnostics require one independent control per design cell")
        self.hydraulics = hydraulics
        self.graph = graph if graph is not None else ConnectivityGraph(hydraulics.pipeline.grid)

    def certificate(self, design):
        return self.graph.evaluate_density(self.hydraulics.pipeline, self.hydraulics.density(design))

    def evaluate(self, design):
        flow = self.hydraulics.evaluate(design)
        certificate, jacobian = self.certificate(design)
        scale = self.graph.minimum_conductance
        return replace(flow, constraints=np.r_[flow.constraints, certificate.constraints / scale],
                       jacobian=np.vstack((flow.jacobian, jacobian / scale)))


def run_connectivity_verification(spacing_mm=0.5, archive=None):
    """Exercise graph cut implications and independently perturb the forward certificate."""
    from .density import DensityPipeline, DesignGrid
    from .reconstruction import reconstruct

    grid = DesignGrid(spacing_mm)
    pipeline = DensityPipeline(grid)
    graph = ConnectivityGraph(grid)
    random = np.random.default_rng(432996)
    samples = {}

    def register(projected):
        answer = graph.evaluate(projected)
        fitted = reconstruct(grid, projected)
        accepted = bool(np.all(answer.constraints <= 0.0))
        connected = fitted.port_component is not None
        return {"conductance": answer.conductance.tolist(), "constraints": answer.constraints.tolist(),
                "hard_paths": list(answer.hard_paths), "certificate_pass": accepted,
                "gamma_half_connected": connected, "port_components": list(fitted.port_components)}

    samples["sharp"] = register(grid.base_raw)
    projected_sharp = pipeline.evaluate(grid.base_raw, q=.01, beta=0, alpha_max=1e5).projected
    samples["filtered_sharp"] = register(projected_sharp)
    barrier = np.full(grid.size, .82)
    col = np.abs(grid.centers_mm[:, 0] - (13.0 + spacing_mm / 2.0)) < spacing_mm / 4.0
    barrier[col] = 0.0
    samples["blocked_main"] = register(barrier)
    gray = barrier.copy()
    gray[col] = .49
    samples["gray_trap"] = register(gray)
    # A leaky auxiliary model with no exact cutoff rates the same failed
    # reconstruction as connected. This never decides certificate feasibility.
    naive_gray_scores, _, _, _ = graph._solve(np.r_[gray, 1., 1., 1.], np.zeros(grid.size + 3))
    samples["gray_trap"]["unguarded_conductance"] = naive_gray_scores.tolist()

    if archive is not None:
        record = json.loads(Path(archive).read_text(encoding="utf-8"))
        if not record["completed"] or record["results"]["algorithm"] != "mma":
            raise ValueError("Expected a completed archived MMA endpoint")
        design = np.asarray(record["results"]["final_checkpoint"]["design"], dtype=float)
        projected = pipeline.evaluate(design, q=.01, beta=0, alpha_max=1e5).projected
        samples["archived_mma"] = register(projected)

    projected = random.uniform(.56, .63, grid.size)
    raw = random.uniform(.56, .63, grid.size)
    direction = random.uniform(-.4, .4, grid.size)
    cotangent = np.array([.9, -.4])
    state = pipeline.evaluate(raw, q=.01, beta=2, alpha_max=1e5)
    answer, raw_jacobian = graph.evaluate_density(pipeline, state)
    filtered_direction = pipeline.filter(direction)
    projected_action = answer.constraint_projected_jacobian @ (state.projection_derivative * filtered_direction)
    dot_error = abs(cotangent @ projected_action - (cotangent @ raw_jacobian) @ direction)

    def raw_constraints(value):
        density = pipeline.evaluate(value, q=.01, beta=2, alpha_max=1e5)
        return graph.evaluate(density.projected).constraints

    projected_action = graph.evaluate(projected).constraint_projected_jacobian @ direction
    raw_action = raw_jacobian @ direction
    fd = []
    for step in (1e-3, 1e-4, 1e-5, 1e-6):
        raw_fd = (raw_constraints(raw + step * direction) - raw_constraints(raw - step * direction)) / (2 * step)
        projected_fd = (graph.evaluate(projected + step * direction).constraints
                        - graph.evaluate(projected - step * direction).constraints) / (2 * step)
        fd.append({"step": step,
                   "raw_inf_error": float(np.linalg.norm(raw_action - raw_fd, np.inf)),
                   "projected_inf_error": float(np.linalg.norm(projected_action - projected_fd, np.inf))})

    derivative_limit = 1e-6
    checks = {"sharp_pass_and_connected": samples["sharp"]["certificate_pass"] and samples["sharp"]["gamma_half_connected"],
              "every_pass_connected": all(not row["certificate_pass"] or row["gamma_half_connected"] for row in samples.values()),
              "solid_cut_rejected": not samples["blocked_main"]["hard_paths"][0]
                                    and samples["blocked_main"]["conductance"][0] <= graph.floor_budget,
              "gray_cut_rejected": bool(not samples["gray_trap"]["gamma_half_connected"]
                                       and not samples["gray_trap"]["certificate_pass"]
                                       and samples["gray_trap"]["conductance"][0] <= graph.floor_budget
                                       and naive_gray_scores[0] >= graph.minimum_conductance),
              "adjoint_dot": bool(dot_error <= 1e-11),
              "directional_fd": sum(max(row["raw_inf_error"], row["projected_inf_error"]) <= derivative_limit
                                    for row in fd) >= 2}
    if archive is not None:
        checks["archived_failure_reproduced"] = (not samples["archived_mma"]["hard_paths"][0]
                                                 and not samples["archived_mma"]["gamma_half_connected"])
    return {"passed": all(checks.values()), "checks": checks, "samples": samples,
            "threshold": graph.threshold, "cutoff": graph.cutoff, "transition_width": graph.transition_width,
            "floor_budget": graph.floor_budget, "minimum_conductance": graph.minimum_conductance,
            "seam_counts_inlet_main_side": graph.seam_counts, "spacing_mm": spacing_mm, "cell_count": grid.size,
            "projected_and_raw_fd": fd, "dot_error": float(dot_error), "job_id": os.environ.get("SLURM_JOB_ID")}


def run_baffle_connectivity_verification(spacing_mm=0.5):
    """Check exact obstacle exclusion, clearance derivatives and fitted connectivity."""
    from shapely.geometry import LineString, box

    from .density import DesignGrid
    from .reconstruction import _PORTS_MM, reconstruct

    grid = DesignGrid(spacing_mm)
    baffle = RectangularBaffleClearance()
    graph = ConnectivityGraph(grid, baffle=baffle)
    sharp = grid.base_raw
    clipped_roof = sharp.copy()
    y = grid.centers_mm[:, 1]
    clipped_roof[(y > 3.5) & (grid.centers_mm[:, 0] > 12.25)] = 0.

    def check_reconstruction(projected, height):
        shape = reconstruct(grid, projected).geometry
        if height > 0.:
            shape = shape.difference(box(baffle.x_min_mm, baffle.wall_y_mm, baffle.x_max_mm,
                                         baffle.wall_y_mm + height))
        polygons = (shape,) if shape.geom_type == "Polygon" else shape.geoms
        polygons = [polygon for polygon in polygons if polygon.geom_type == "Polygon"]
        ports = [tuple(i for i, polygon in enumerate(polygons) if polygon.boundary.covers(LineString(port)))
                 for port in _PORTS_MM]
        return [bool(ports[0] and ports[0] == port) for port in ports[1:]]

    # Use an actual row's clearance-transition midpoint. A horizontal graph
    # edge on that row crosses x=13 or x=14 and is neither laterally nor
    # base-clear. This makes its direct and effective derivatives nontrivial.
    row = float(np.unique(y[y >= .5])[0])
    center_height = row - baffle.margin_mm - .5 * baffle.transition_mm
    edge = np.array([[[12.75, row], [13.25, row]]])
    expected_gate = .5
    expected_edge_derivative = -1.875 / baffle.transition_mm
    edge_value, edge_derivative = baffle.evaluate(edge, center_height)
    physical_heights = (0., center_height, 1.5, baffle.max_height_mm)
    cases = []
    for label, projected in (("sharp", sharp), ("roof_restricted", clipped_roof)):
        for height in physical_heights:
            answer = graph.evaluate(projected, height)
            clipped_connected = check_reconstruction(projected, height)
            cases.append(dict(fixture=label, height_mm=height, conductance=answer.conductance.tolist(),
                              hard_paths=list(answer.hard_paths), post_subtraction_paths=clipped_connected,
                              conductance_pass=(answer.constraints <= 0.).tolist()))

    def conductance(height):
        return graph.evaluate(sharp, height).conductance

    checks = []
    for label, height, mode in (("lower_bound", 0., "forward"), ("transition", center_height, "centered"),
                                ("upper_bound", baffle.max_height_mm, "backward")):
        result = graph.evaluate(sharp, height)
        fd_report = []
        for step in (1e-3, 1e-4, 1e-5, 1e-6):
            if mode == "centered":
                fd = (conductance(height + step) - conductance(height - step)) / (2 * step)
            else:
                sign = 1. if mode == "forward" else -1.
                fd = sign * (-3. * result.conductance + 4. * conductance(height + sign * step)
                             - conductance(height + sign * 2. * step)) / (2 * step)
            fd_report.append(dict(step_mm=step, derivative=fd.tolist(),
                                  absolute_error=float(np.max(np.abs(fd - result.height_jacobian)))))
        checks.append(dict(label=label, height_mm=height, scheme=mode, height_jacobian=result.height_jacobian.tolist(),
                           differences=fd_report,
                           passed=sum(f["absolute_error"] <= 1e-6 for f in fd_report) >= 2))

    case_checks = dict(no_obstructed_edge_pass=all(all(not certified or actual
                                                       for certified, actual in zip(case["hard_paths"],
                                                                                    case["post_subtraction_paths"]))
                                                   for case in cases),
                       sharp_half_height_leaves_both_outlets_connected=all(cases[2]["conductance_pass"]),
                       high_obstruction_rejected=not cases[-1]["post_subtraction_paths"][0]
                                                   and not cases[-1]["conductance_pass"][0],
                       gate_transition_value=bool(abs(edge_value[0] - expected_gate) <= 1e-13),
                       gate_transition_derivative=bool(abs(edge_derivative[0] - expected_edge_derivative) <= 1e-12),
                       nonzero_effective_transition_derivative=max(map(abs, checks[1]["height_jacobian"])) > 1e-8)
    return dict(passed=all(case_checks.values()) and all(check["passed"] for check in checks),
                spacing_mm=spacing_mm, graph_checks=case_checks, fixtures=cases,
                height_derivative=checks, job_id=os.environ.get("SLURM_JOB_ID"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--spacing", type=float, default=.5)
    parser.add_argument("--archive", type=Path, help="Optional archived 240-cell MMA summary")
    parser.add_argument("--baffle", action="store_true", help="Verify the additional rectangular baffle gate")
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this numerical verification as a Sol Slurm job")
    report = run_connectivity_verification(args.spacing, args.archive)
    if args.baffle:
        report["baffle"] = run_baffle_connectivity_verification(args.spacing)
        report["passed"] = report["passed"] and report["baffle"]["passed"]
    with args.output.open("x", encoding="utf-8") as output:
        json.dump(report, output, indent=2, allow_nan=False)
        output.write("\n")
    raise SystemExit(0 if report["passed"] else 1)


if __name__ == "__main__":
    main()
