# 2D BFS2d_1600 Navier–Stokes recovery test

`problems_navierStokes_2D_3D_bfs_recovery` runs on four MPI ranks. It reuses the
inputs and executable from `../navierStokes_2D_3D_bfs` and generates its reference
solutions independently in an isolated working directory.

This test exercises these operations:

1. Run the uninterrupted case and write its ordinary reference checkpoints.
2. Restart at 0.005 with recovery enabled. A zero linear iteration fraction
   deliberately triggers protective checkpointing for converged real solves.
   Three permitted Newton updates also exercise convergence on the final update.
   Continue to 0.01 and compare the final solution to the reference.
3. Restart from both the latest recovery generation and the protected generation
   preceding the warning, and compare against the same reference at 0.01.
4. Restart at 0.005 with only one Newton iteration permitted and cancellation
   enabled. Require a controlled stop and a complete recovery checkpoint.
5. Restart from that checkpoint with the normal solver settings and compare
   against the reference.
6. Repeat the one-iteration case with cancellation disabled. Require continuation
   to 0.01 while the recovery state stays at 0.005. Restart from this preserved
   state and compare against the reference.

The zero fraction is a test setting; production defaults to `maxIter/2`. The
reference comparisons retain the existing restart error tolerance. The tests
also check completion markers, recovery diagnostics and retention of exactly
one protected and one latest generation.

`problems_recoveryCheckpoint_MPI_2` in `UnitTest` checks the exact threshold boundary,
zero linear solves, final-iteration convergence classification, nonfinite and
failed-linear-solve classification, independent vector copies, and an interrupted
checkpoint write on two MPI ranks. A partial write must leave the previous
complete generation selected by `Latest.xml`.
