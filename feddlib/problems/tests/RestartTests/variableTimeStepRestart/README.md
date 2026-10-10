# Restart with a changed timestep

Four four-rank cases reuse the 2D BFS, Turek FSI, nonlinear-elasticity Newmark
and 3D arterial segment (Absorbing Paper) drivers. Each checkpoints at 0.005
with dt=0.0025, then
restarts with dt=0.0015 and no original interval schedule. An independent
uninterrupted reference changes its dt at 0.005 through `Timestepping Intervalls`.
The comparison at 0.0095 uses the existing 1e-12 relative error bounds.

The resumed simulation writes a new checkpoint at 0.008 and restarts from that
checkpoint too. Neither restart time lies on the global grid of the new dt.
FSI and Newmark retain their existing checkpoint/update order and advance an
extra reference step to publish the comparison state. Newmark compares runtime
displacement, velocity and acceleration; FSI compares its coupled primary fields.
The artery also compares its pressure boundary history, including the transition
between the source checkpoint and the first completed resumed step.

The `timeStepping` unit test independently checks BDF2 against a quadratic,
interval-boundary landing, a shortened final step, saved incoming dt and the
step counter. Metadata tests retain incompatible-field/history checks and
exercise loading version-1 uniform checkpoints with a changed next timestep,
version-2 nonuniform history and invalid or incomplete saved clocks.

The Navier–Stokes and Turek cases also resume the same ParaView files across
both restarts. The Turek case writes each frame's moving mesh. The independent
reference has visualization disabled to preserve the source series. XMF checks
verify increasing timestamps and a single frame per simulation timestep.
