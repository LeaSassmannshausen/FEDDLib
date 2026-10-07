# Checkpointing and restart

Time-dependent problems can write HDF5 checkpoints and resume from a saved
solution. This branch provides the shared solution/BDF/Newmark checkpoint code
and FSI restart support. SCI/FSCI solvers and material-history extensions are
separate features.

Configure an uninterrupted run with these entries in `Timestepping Parameter`:

```xml
<Parameter name="Checkpointing" type="bool" value="true"/>
<Parameter name="Checkpoint directory" type="string" value="restartCheckpoints"/>
<Parameter name="Number Checkpoints" type="int" value="2"/>
<ParameterList name="Checkpoints">
    <Parameter name="1" type="double" value="0.01"/>
    <Parameter name="2" type="double" value="0.02"/>
</ParameterList>
```

Resume using the same mesh, discretization, variable names, time integration
scheme, and time step size, with a later `Final time`:

```xml
<Parameter name="Checkpointing" type="bool" value="false"/>
<Parameter name="Restart" type="bool" value="true"/>
<Parameter name="Restart directory" type="string" value="restartCheckpoints"/>
<Parameter name="Time step" type="double" value="0.01"/>
```

`Time step` is the physical restart time; `dt` is the time increment. Both
directory parameters default to the working directory. Checkpointing and restart
are disabled by default. Use separate output and input directories if a resumed
run also writes checkpoints, since creating an output file truncates it.

Checkpoint times should coincide with time steps. For BDF/FSI, checkpoints are
written when the next step begins; let the uninterrupted run advance one step
past a checkpoint you need. A multistep restart needs the saved previous steps
as well as the solution at the restart time. The tests restart after several
steps, so that this history is available.

The checkpoint files include `Solution<variable>.h5`, Newmark displacement,
velocity and acceleration, and, for FSI, moving-mesh mass products and geometry.
Dataset names are physical times formatted by `std::to_string`. Transfer the
complete checkpoint directory, rather than just the displacement or velocity
file. The FSI test covers explicit geometry with the Turek benchmark; it does
not validate restart of stateful outlet boundary conditions.

## Validation

With tests enabled and the required Trilinos solvers available:

```sh
cmake --build build --target problems_unsteadyNavierStokes_restart problems_fsi_restart
ctest --test-dir build --output-on-failure -R 'problems_(unsteadyNavierStokes_restart|fsi_restart)'
```

- `unsteadyNavierStokes_restart` runs the 2D `BFS2d_1600.mesh` and 3D
  `BFS3dCC.mesh` cases with MPI rank pairs `4 -> 4` and `4 -> 6`. Each test runs
  an uninterrupted reference on four ranks, writes checkpoints at `0.01` and
  `0.02`, then restarts on four or six ranks at `0.01`
  and compares velocity and pressure at `0.02` with relative tolerance `1e-12`.
  Absolute and relative l2 errors are printed for both fields. The shared
  linear and nonlinear solver tolerances are `1e-12` to resolve differences
  below the comparison bound when the MPI partition changes. The problem
  settings are in `parametersProblem.xml`, linear solver settings in
  `parametersSolver.xml`, the mesh/dimension overrides in
  `parametersProblem2D.xml` and `parametersProblem3D.xml`, and the second-phase
  overrides in `parametersProblem_restart.xml`. Test names include the rank
  pair (`2D_4_to_4`, `2D_4_to_6`, `3D_4_to_4`, `3D_4_to_6`). Each combination has
  its own working directory and generates its own checkpoints. The runner
  prepares the inputs quietly and streams both simulations' normal output;
  use `ctest -V` to display the solver output and restart errors.
- `fsi_restart` runs an uninterrupted simulation and a restarted simulation on
  four MPI ranks, comparing fluid velocity, pressure, and solid displacement at
  the same final time with relative tolerance `1e-12`. Both phases load the case
  from `parametersProblemFSI.xml`; `parametersProblemFSI_restart.xml` contains
  only the overrides to resume at `0.01` and stop at `0.02`. The test disables
  visualization and benchmark exports and reports the relative error for each
  of the three fields.

The restart tests use upstream meshes and do not require the new SCI meshes or
the AceGen Interface2 material models.

## Refactoring in increments

Each increment should keep the uninterrupted/restarted comparisons passing and
retain a negative check in which only an input history value is perturbed.
Do not combine changes to restoration order, numerical history updates, and
the on-disk format in one increment.

| Increment | Deliverable | Acceptance check |
| --- | --- | --- |
| 1a: separate checkpoint reads | Named restore methods for primary fields, BDF history, Newmark displacement/derivatives, FSI mass products, and ALE geometry. This increment retains the existing call order and file format. | Both Navier-Stokes dimensions/rank counts and FSI pass; perturbed BDF history fails. |
| 1b: explicit restore phase | A coordinator invokes component restore methods in a documented order before time advancement. Remove lazy imports from numerical update routines after handling the first-update history shift explicitly. | Preserve first-step and final-step equivalence; exercise both standalone and FSI Newmark timing conventions. |
| 2: accepted-step checkpoints | Capture a consistent state after an accepted timestep, including its complete integration history. | Restart from the final timestep without advancing one extra step; verify BDF, Newmark, and moving geometry. |
| 3: metadata and compatibility | Add a versioned manifest with time, step number, integration settings, mesh/DOF identity, field names/sizes, and required history. Validate before changing simulation state. | Reject wrong mesh/discretization/method and missing or inconsistent history with clear diagnostics. |
| 4: publish complete checkpoints | Write each checkpoint to a new temporary directory, publish only after all MPI ranks/components succeed, and retain the preceding complete checkpoint. | Interrupted writes leave the preceding checkpoint usable; resuming and writing cannot truncate the input checkpoint. |
| 5: timestep identifiers | Use integer step identifiers and history indices as keys; store physical times as full-precision metadata. Define how legacy checkpoints are recognized/read. | Very small dt and non-integer time values do not collide; supported legacy checkpoints remain readable. |

Increment 1a separates file reads from the surrounding numerical operations:

- `Problem::initializeVectors()` allocates vectors only. Both linear and
  nonlinear `initializeProblem()` retain automatic restoration through
  `restoreSolutionFromCheckpoint()` when restart is enabled. Explicit callers
  must restore before assembly or linking subproblem vectors.
- `TimeProblem` has separate restore methods for multistep history, moving-mesh
  mass products, Newmark displacement history, and Newmark derivatives.
- `FSI` restores its old/new ALE displacement buffers through a named method.

Restoration is not yet a single coordinated phase. The legacy Newmark format
labels different states depending on whether the caller's clock denotes the
start or end of a step. The first subsequent update must also avoid shifting
restored BDF history or advancing Newmark derivatives twice. Those dependencies
are the scope of increment 1b; existing timing and update formulas are retained
in 1a.

Validation of increment 1a: all five restart cases (2D/3D Navier-Stokes on four
and six ranks, plus FSI on four ranks) passed at relative tolerance `1e-12`,
with zero reported comparison error. In a disposable copy of the 2D four-rank
checkpoint, scaling only velocity history at `0.008` by `1.1` left the fields
at the restart time and the final reference unchanged. The restarted run
failed as expected, with relative velocity error approximately `0.0103`.
