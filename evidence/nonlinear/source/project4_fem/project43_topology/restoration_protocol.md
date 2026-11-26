# Gate 2 optimizer recovery

The existing mathematical review conversation accepted the tested density, discrete adjoint and reconstruction implementations. It authorized these restoration experiments while leaving production approval pending.

Phase I uses SLSQP to minimize the nonnegative epigraph slack \(t\). The constraints are \(g_i(m)+\delta-t\leq0\) with the design in the unit box. Initialize \(t_0=1+\max(0,\max_i[g_i(m_0)+\delta])\), and use the original binary topology or all-zero reduced coefficients. The physical design evaluation is cached by \(m\), even if the slack changes.

Start with a scaled margin of 0.05, which targets 49.905–50.095% side flow and at most 82.95 mm². Use up to 250 unique physical evaluations per attempt. If the tighter restoration attempt fails, preserve it and try the original constraints with zero margin. Verify the engineering and margin feasibility independently before handing off an evaluated design. Archive the optimizer endpoint even when a different evaluated trial supplies the feasible start.

Phase II uses the original pumping-power objective and constraints. In the 240-cell verification fixture, run SLSQP and MMA/CCSA from one common restored design with identical additional 250-evaluation budgets. Retain the initial binary-design pilots as separate evidence. Independently verify the terminal designs using actual projected KKT residuals. Separately attempt restoration in the 20-variable space. If feasible, transfer that design to the full space to verify matching responses.

Keep the original \(P_0\) throughout continuation. Reevaluate feasibility after a stage changes and repeat restoration if needed. Count restoration, forward solves, adjoint right-hand sides and Newton iterations in cost comparisons. Process peak RSS is a high-water mark, so report it with that scope.

The independent KKT witness minimizes a bounded stationarity slack subject to the same stationarity and complementarity inequalities. It recomputes the projected residual at the original design coordinates. The review kept all four KKT thresholds and required independent analytic fixtures for active, near-active and interior variables.
