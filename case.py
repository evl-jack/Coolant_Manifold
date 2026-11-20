"""Editable inputs for the first 2D RANS verification case (SI units)."""

from mesh import BRANCH_ROOT, FIXED_MAIN_RUN, FIXED_SIDE_RUN, MAIN_CHANNEL, SIDE_BRANCH

CASE = {
    "name": "straight_channel_baseline",
    "main_channel": MAIN_CHANNEL,
    "side_branch": SIDE_BRANCH,
    "branch_root": BRANCH_ROOT,
    "fixed_main_run_m": FIXED_MAIN_RUN,
    "fixed_side_run_m": FIXED_SIDE_RUN,
    "include_side_outlet": False,        # first solver baseline is a rectangle
    "spacing_m": 0.001,
    "density_kg_m3": 998.2,             # liquid water, nominal 20 C
    "dynamic_viscosity_pa_s": 1.002e-3, # liquid water, nominal 20 C
    "inlet_velocity_m_s": (2.0, 0.0),
    "outlet_pressure_pa": 0.0,          # gauge pressure at both outlets
    "side_to_main_ratio": 0.75,
    "design_volume_fraction": 0.30,    # fluid-equivalent fraction of editable area
    "outlet_balance_tolerance": 0.001, # allowable mismatch as a fraction of inlet flow
    "turbulence_intensity": 0.05,       # starter assumption; replace with measured value
    "turbulence_length_scale_m": 0.00056, # 0.07 x 8 mm hydraulic diameter
    "max_iterations": 500,
    "velocity_relaxation": 0.5,
    "pressure_relaxation": 0.3,
    "velocity_correction_relaxation": 0.5,
    "convergence_tolerance": 1e-5,
    "pseudo_cfl": 2.0,
}
