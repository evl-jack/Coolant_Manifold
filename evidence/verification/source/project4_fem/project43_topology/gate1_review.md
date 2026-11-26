# Gate 1 mathematical review

The source audit and implementation plan were submitted through Chrome to the existing [Brinkman SST Modeling conversation](https://chatgpt.com/g/g-p-5df5c39440c08191bf518d6456de7dbc-jacks-solver/c/6b05bfb1-9b00-83e8-8d32-af9075dbf801). The reviewer approved implementation and Gate 2 verification. Production optimization requires separate approval after Gates 2 and 3.

## Decisions adopted

- Use planar Taylor–Hood FEM and damped Newton on the exact 101 mm² computational union. Fix the geometry and tags of the three external ports. Use the existing signed port-power calculation in W/m and a shared fixed baseline reference.
- Evaluate resistance on the projected density of each design cell. Use piecewise-constant transfer on a conforming analysis mesh and zero resistance in fixed fluid. Check the transfer transpose to a relative error of 1e-12. Propagate the adjoint with constrained DOFs set to zero and include quadrature weights, the resistance scale, and the filter transpose.
- Center-sample the Wendland controls and retain group normalization. Report singular values, conditioning, column norms, and maximum changes per cell, with spatial expansion and contraction coverage. Test reduced-to-full equivalence for downstream quantities and for gradients pulled back with the reduced map transpose.
- Keep homogeneous Neumann filter boundaries. Compare the sharp fitted baseline, unfiltered binary immersed baseline, and filtered binary-initialized immersed baseline. Record seam density and resistance at x = 6 mm, x = 16 mm, and along the 1 mm-wide side seam at y = 6 mm.
- Use SLSQP and NLopt LD_MMA, identified as MMA/CCSA. Verify both on a constrained problem with a known solution and compare both on the 240-variable topology pilot from the prescribed infeasible start. Verify the installed version and constraint signs on Sol.
- For the three-multiplier KKT LP, enforce projected stationarity directly at the actual box coordinates, without snapping. Recompute feasibility, stationarity, complementarity and dual feasibility from the recovered witness.
- Freeze an unsmoothed 0.5 contour from bilinear interpolation with constant-normal extension to the design boundary. Use a bilinear asymptotic decider for marching squares, treat degenerate point-only contacts as disconnected, and union with exact fixed fluid. Archive both the contour and polygon. Include isolated pockets in geometric area, give each solved closed pocket its own pressure gauge, and retain solid islands as no-slip boundaries.
- Use at least three analysis levels in the final resistance screen and fitted-wall assessment. Examine both successive changes. Certify physical flow feasibility only if its entire uncertainty interval lies in [0.499, 0.501]. Certify area only if the upper uncertainty bound is at most 83 mm².
- Keep cached evaluations scoped to the mesh and stage, and record failed Newton trials without presenting an unconverged result as objective or gradient data.

## Gate 3 limits for unfiltered binary fixtures

| Measure | Maximum |
|---|---:|
| Gross solid ingress / inlet flow | 0.10% |
| Gross solid egress / inlet flow | 0.10% |
| Absolute through-solid flux on each predefined section / inlet flow | 0.02% |
| Solid RMS speed / bulk inlet speed | 1.0% |
| Immersed versus fitted flow split | 0.02 percentage points |
| Immersed versus fitted power | 0.5% |

For filtered cases, report the leakage diagnostics for the 0.5 nominal solid region and the density-at-most-0.05 solid core. Include intermediate density. Gate 3 remains pending until a verified binary mesh and resistance combination satisfies the prospective criteria.
