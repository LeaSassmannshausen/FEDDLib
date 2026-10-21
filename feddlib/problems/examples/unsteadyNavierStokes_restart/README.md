# Unsteady Navier-Stokes with restart

This example derives from `../unsteadyNavierStokes` and keeps only the 3D
unstructured mesh path. The supplied `BFS3dCC.mesh` is a backward-facing step.
Velocity uses P2 elements and pressure uses P1; P1/P1 is also supported.
Boundary markers are 1 for walls, 2 for the prescribed parabolic inlet and 3 for
natural outflow; marker 4, when present, is treated as a no-slip obstacle.
Replacing the mesh requires matching these boundary markers and the inlet profile.

`Parameter/Inflow ramp time` sets a smooth half-cosine startup ramp in seconds
(default: 0.01 s). The inlet starts at zero, reaches half its final velocity at
half the ramp time, then reaches `MaxVelocity` at the ramp time and stays steady.
Set the ramp time to zero to apply the full profile immediately. The ramp uses
physical simulation time, so a restart during the ramp continues from its saved
time. Keep the ramp setting unchanged between the producer and restart.

Build the target from the workspace directory containing `source/` and `build/`:

```sh
cmake --build build --target problems_unsteadyNavierStokes_restart -j 2
cd build/feddlib/problems/examples/unsteadyNavierStokes_restart
```

Run the producer from 0 to 0.01 s, then continue from its final checkpoint to 0.02 s:

```sh
mpirun -np 4 ./problems_unsteadyNavierStokes_restart.exe
mpirun -np 4 ./problems_unsteadyNavierStokes_restart.exe --restartfile=parametersProblem_restart.xml
```

The base problem, preconditioner and solver files default to `parametersProblem.xml`,
`parametersPrec.xml` and `parametersSolver.xml`. Override their paths with
`--problemfile`, `--precfile` and `--solverfile`. `--restartfile` merges an optional
override into the base problem settings. Relative input/output paths are resolved
from the launch working directory.

The producer saves `Solutionu.h5`, `Solutionp.h5` and checkpoint metadata in
`checkpoints/`, including the previous solution needed by BDF2. `Time step` selects
the physical time to restore; `dt` is the step duration. Restoration validates the
mesh, DOFs, fields and integration history before continuing. Navier-Stokes writes
checkpoints after accepted steps, so the final checkpoint needs no extra timestep.
The continuation writes its final checkpoint to `continuedCheckpoints/`; never
use the restart input directory as the destination for new checkpoints.

`Exporter/Resume output=true` in the restart override continues the producer's
ParaView (`u` and `p`) and diagnostic output in the same working directory. To start
a new output series instead, set it to false; `Keep old output=true` archives the
previous files under `old_output/run_N`. The fresh configuration already uses
this archive option. View the `.xmf` files in ParaView.

The same unstructured input mesh can be repartitioned on a different MPI rank
count for restart. Keep its global geometry, numbering and discretization
compatible with the checkpoint. In other drivers using **structured partitions**,
checkpointing now warns that changing MPI rank count or H/h can change the mesh
and DOFs. Use the same rank count, H/h and geometry settings to reproduce those
meshes. This example rejects structured meshes and dimensions other than 3.

Solver timings and checkpoint feature timings are reported at the end, including
the stacked timer report. Disable `General/Checkpoint timings` to omit the
checkpoint feature counters.
