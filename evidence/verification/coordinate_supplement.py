from pathlib import Path
import project4_fem.project43_topology.verification.fixed_pipe_adjoint_audit as a
a.STEP_SIZES = (.04, .02, .01, .005, .002, .001, .0001, .00001, .000001)
result = a.audit(Path("input"), ("beta_0", "beta_2", "beta_8"), ("junction_cell", "collar_cell"), Path("adjoint-coordinate-supplement.json"))
raise SystemExit(not result["passed"])
