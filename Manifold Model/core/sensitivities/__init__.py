"""Local derivatives and residual-based sensitivity building blocks.

These helpers do not supply missing derivatives of a forward CFD model.
No total flow, geometry, or RANS gradient is claimed until it is assembled
and checked with that model.
"""
