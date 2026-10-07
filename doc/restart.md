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

Checkpoint times should coincide with time steps. Standalone BDF problems,
including Navier-Stokes, write checkpoints after a successful solve and clock
advancement. A run ending at a requested checkpoint time writes that checkpoint
without an extra step. The current solution and required previous solutions are
saved without shifting numerical history a second time. This covers fixed-step
BDF integration and retains the existing solution files and time keys.

FSI and standalone Newmark retain their existing checkpoint timing. For both,
let the uninterrupted run advance one step past a checkpoint you need. A
multistep restart needs the saved previous steps as well as the solution at the
restart time. The tests restart after several steps so this history is available.

The checkpoint files include `Solution<variable>.h5`, Newmark displacement,
velocity and acceleration, and, for FSI, moving-mesh mass products and geometry.
Dataset names are physical times formatted by `std::to_string`. Transfer the
complete checkpoint directory, rather than just the displacement or velocity
file. The FSI test covers explicit geometry with the Turek benchmark; it does
not validate restart of stateful outlet boundary conditions.

## Validation

With tests enabled and the required Trilinos solvers available:

```sh
cmake --build build --target problems_unsteadyNavierStokes_restart problems_fsi_restart problems_unsteadyNonLinElasticity_restart
ctest --test-dir build --output-on-failure -R 'problems_(unsteadyNavierStokes_restart|fsi_restart|unsteadyNonLinElasticity_restart)'
```

- `unsteadyNavierStokes_restart` runs the 2D `BFS2d_1600.mesh` and 3D
  `BFS3dCC.mesh` cases with MPI rank pairs `4 -> 4` and `4 -> 6`. Each test runs
  an uninterrupted reference on four ranks, writes checkpoints at `0.01` and
  `0.02`, stops exactly at `0.02`, then restarts on four or six ranks at `0.01`
  and compares velocity and pressure at `0.02` with relative tolerance `1e-12`.
  Absolute and relative l2 errors are printed for both fields. The shared
  linear and nonlinear solver tolerances are `1e-12` to resolve differences
  below the comparison bound when the MPI partition changes. The problem
  settings are in `parametersProblem.xml`, linear solver settings in
  `parametersSolver.xml`, the mesh/dimension overrides in
  `parametersProblem2D.xml` and `parametersProblem3D.xml`, and the second-phase
  overrides in `parametersProblem_restart.xml`. Each test also runs an
  independent reference stopping at `0.03` (`parametersProblem_reference.xml`),
  then restarts the first run's final checkpoint at `0.02` and continues to
  `0.03` (`parametersProblem_restart_final.xml`). This comparison reads from
  `referenceCheckpoints/`, while restoration reads from `restartCheckpoints/`,
  so the extended reference cannot supply a missing producer checkpoint.
  Both comparisons use the same error bound. Test names include the rank
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

- `unsteadyNonLinElasticity_restart` uses the 2D P2 `square_solid.mesh` case
  from the `unsteadyNonLinElasticity` example, with Saint Venant-Kirchhoff
  material, a clamped left edge and a constant surface load on the right edge.
  Newmark uses `beta = 0.25`, `gamma = 0.5` and `dt = 0.0025`. The uninterrupted
  run writes checkpoints at `0.01` and `0.02`; the second phase resumes at
  `0.01`. Rank pairs are `4 -> 4` and `4 -> 6`, in separate working directories.
  Both phases stop at `0.0225` because the current Newmark loop finalizes and
  writes the `0.02` state at the beginning of the following step. The test
  compares runtime displacement history, velocity and acceleration at `0.02`
  with the reference checkpoint, prints absolute/relative l2 errors and requires
  relative errors at most `1e-12`. Each reference field must be nonzero.
  Newmark checkpoint timing is deliberately retained for this baseline test.
  Standalone nonlinear Newmark now also writes the primary displacement file
  required by initialization; elasticity assembly preserves a restored field.
  Both rank pairs pass at `1e-12`; the largest relative error in the initial
  validation was approximately `4.66e-14` (acceleration, `4 -> 6`). In separate
  disposable checkpoints, scaling only velocity or only acceleration at `0.01`
  by `1.1` makes continuation fail the numerical comparison. The FSI regression
  also passes with these changes.

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

Increment 2 for standalone BDF: the linear and nonlinear multistep loops call
`TimeProblem::writeMultistepCheckpoint(completedTime, dt)` after a successful
step. Their history updates disable the legacy start-of-step export. The writer
stores the current solution first, followed by the required unshifted history;
it includes time-zero history for a first-step BDF2 restart and writes shared
history entries only once for adjacent checkpoints. It does not change the
solution, history or clock. `Safe all solution` retains
initial-state output and also includes the final accepted solution. FSI and
Newmark still use the legacy path and need a separate increment to finalize
derivatives and moving-mesh state before writing a coupled checkpoint.

Validation of the BDF increment: all four 2D/3D rank-pair cases pass both the
intermediate and final-checkpoint continuation comparisons at `1e-12`. The
largest relative error was approximately `1.60e-13`. The FSI regression and
isolated first-step/adjacent-checkpoint and save-all checks also pass. Scaling
only velocity history at `0.0175` by `1.1`, while leaving the restart state at
`0.02` and the independent reference unchanged, makes continuation fail with
relative velocity error approximately `0.0194`.
