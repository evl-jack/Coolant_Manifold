"""Time five matched fixed-pipe MMA evaluations on Windows."""

import ctypes
from ctypes import wintypes
import json
from pathlib import Path
import statistics
import sys
from time import perf_counter

import numpy as np

from project4_fem.project43_topology.evaluator import Stage
from project4_fem.project43_topology.fixed_pipe_pilot_model import FixedPipeEvaluator, build_fixed_pipe
from project4_fem.project43_topology.gate2 import environment


class ProcessMemory(ctypes.Structure):
    _fields_ = [("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD)] + [
        (name, ctypes.c_size_t) for name in (
            "PeakWorkingSetSize", "WorkingSetSize", "QuotaPeakPagedPoolUsage", "QuotaPagedPoolUsage",
            "QuotaPeakNonPagedPoolUsage", "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage",
        )
    ]


def process_memory():
    current = ctypes.windll.kernel32.GetCurrentProcess
    current.restype = wintypes.HANDLE
    memory = ProcessMemory()
    memory.cb = ctypes.sizeof(memory)
    get_memory = ctypes.windll.psapi.GetProcessMemoryInfo
    get_memory.argtypes = [wintypes.HANDLE, ctypes.POINTER(ProcessMemory), wintypes.DWORD]
    if not get_memory(current(), ctypes.byref(memory), memory.cb):
        raise ctypes.WinError()
    return {"rss_MiB": memory.WorkingSetSize / 2**20, "peak_rss_MiB": memory.PeakWorkingSetSize / 2**20}


def main():
    fixture = json.loads(Path(sys.argv[1]).read_text(encoding="utf-8"))
    seeds = np.asarray(fixture["designs"], dtype=float)
    domain, pipeline = build_fixed_pipe()
    stage = Stage(.10, 0., 1e5)
    reference = FixedPipeEvaluator(domain, pipeline, stage, 1., 0.)
    baseline = reference.evaluate(seeds[0])
    ev = FixedPipeEvaluator(domain, pipeline, stage, baseline.power_w_per_m, baseline.side_fraction,
                            initial_state=reference.state)
    ev.evaluate(seeds[0])
    rows = []
    for number, design in enumerate(seeds[1:], start=2):
        forward, adjoint = ev.forward_seconds, ev.adjoint_seconds
        start = perf_counter()
        result = ev.evaluate(design)
        elapsed = perf_counter() - start
        row = {"sol_mma_unique_evaluation": number, "elapsed_seconds": elapsed,
               "forward_seconds": ev.forward_seconds - forward, "adjoint_seconds": ev.adjoint_seconds - adjoint,
               "normalized_power": result.objective, "side_fraction": result.side_fraction,
               "nonlinear_residual": result.nonlinear_residual, **process_memory()}
        rows.append(row)
        print(json.dumps(row), flush=True)
    if ev.forward_solves != 6 or ev.adjoint_rhs != 12:
        raise RuntimeError("Expected five timed and one warm-up forward solutions and two adjoints per call")
    average = statistics.mean(row["elapsed_seconds"] for row in rows)
    sol_mean = statistics.mean(fixture["sol_evaluation_seconds"])
    summary = {"environment": environment(), "stage": {"q": .10, "beta": 0., "alpha_max": 1e5},
               "evaluations": rows, "average_seconds": average, "peak_rss_MiB": max(r["peak_rss_MiB"] for r in rows),
               "fem_triangles": len(domain.fem.triangles), "fem_dofs": domain.fem.size,
               "free_fem_dofs": len(ev.state.free), "editable_variables": pipeline.size,
               "sol_job": fixture["sol_job"], "sol_average_seconds": sol_mean,
               "local_over_sol_ratio": average / sol_mean, "estimated_300_evaluation_seconds": 300 * average,
               "initial_normalized_power": 1., "initial_side_fraction": baseline.side_fraction}
    Path(sys.argv[2]).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k != "evaluations"}), flush=True)


if __name__ == "__main__":
    main()
