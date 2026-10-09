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

New checkpoints include a version-1 XML manifest for each component and time,
for example `Checkpoint_u_p_0.010000.xml`. Component names distinguish a coupled
FSI problem from its fluid and structure subproblems. The manifest records the
format/layout, physical time, step number, fixed `dt`, BDF order and extrapolation
history, Newmark parameters where applicable, field names, FE types, components,
global DOF counts/index bases, and the exact required HDF5 files and time keys.

`core/Checkpointing/CheckpointMetadata.hpp` owns schema construction (`makeSchema`),
metadata keys, integration/history compatibility rules, and manifest I/O.
`core/Checkpointing/CheckpointMeshFingerprint.hpp` computes the reference mesh/DOF identity.
`Problem::prepareCheckpointMetadata()` collects settings and field descriptions,
caches the schema before ALE motion, and triggers restart validation.

Before primary fields, geometry or history are restored, the reader compares
this metadata and inspects every required dataset. Mesh/DOF identity is checked
using two order-independent 64-bit record checksums over global node IDs,
reference coordinates, node markers, element IDs/connectivity/markers,
boundary subelements and global field DOF IDs. These are compatibility
fingerprints, not cryptographic integrity checks. They ignore MPI ownership,
so rank count may change while global numbering must remain consistent.
The reference fingerprint is cached before ALE moves the mesh. Solver
tolerances, output settings and final time may change.

Missing manifests are rejected by default. To explicitly read an older
checkpoint, set `Allow legacy restart` to `true` in `Timestepping Parameter`.
Legacy mode still checks the required datasets and dimensions; it cannot verify
mesh or integration compatibility. Version 1 supports fixed-step BDF1/BDF2,
Newmark, and the existing FSI BDF/Newmark layouts. Segmented timestep schedules
are rejected when checkpointing/restart is enabled. Numerical history updates
and the Newmark/FSI checkpoint timing remain unchanged.

The XML manifest is descriptive metadata, **not** a completion marker. The
reader checks the actual HDF5 data even when a manifest exists. Publishing whole
checkpoints atomically and replacing decimal dataset keys remain separate work.

`problems_checkpointMetadata_MPI_2` exercises incompatible mesh/FE/DOF metadata,
BDF order and timestep changes, unknown versions, missing history, inconsistent
HDF5 length/vector-count/shape/type, a mismatched destination map, and explicit
legacy mode. Rejections are checked on every MPI rank. The 2D 4-to-4
Navier–Stokes restart test additionally rejects a same-size mesh with one changed
coordinate and modified FE/BDF metadata through the actual restore entry point.
The metadata unit test also checks initial solution compatibility, including new
integration settings, absent history, required primary fields and rejection of
legacy manifests, FSI and Newmark initialization.

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
FSI copies the coupled `Timestepping Parameter` settings into its fluid and solid
components before constructing them. Configure checkpointing/restart on the
coupled problem; its dt, checkpoint schedule and directories also govern fluid
mass history and Newmark displacement, velocity and acceleration. Driver code
does not need to copy these settings into the component parameter lists.

The checkpoint files include `Solution<variable>.h5`, Newmark displacement,
velocity and acceleration, and, for FSI, moving-mesh mass products and geometry.
Dataset names are physical times formatted by `std::to_string`. Transfer the
complete checkpoint directory, rather than just the displacement or velocity
file. The FSI test covers explicit geometry with the Turek benchmark; it does
not validate restart of stateful outlet boundary conditions.

## Initial solution for a new Navier–Stokes simulation

Standalone multistep Navier–Stokes can load saved velocity and pressure as the
initial condition of a new simulation. Configure `Timestepping Parameter` with:

```xml
<Parameter name="Restart" type="bool" value="false"/>
<Parameter name="Initial solution" type="bool" value="true"/>
<Parameter name="Initial solution directory" type="string" value="initialSource"/>
<Parameter name="Initial solution time" type="double" value="0.015"/>
<Parameter name="Class" type="string" value="Multistep"/>
<Parameter name="BDF" type="int" value="2"/>
<Parameter name="dt" type="double" value="0.004"/>
<Parameter name="Final time" type="double" value="0.008"/>
```

`Initial solution` defaults to false. When enabled, its directory and source time
must be specified. The source time identifies `Solutionu.h5`, `Solutionp.h5` and
their versioned manifest; it does not set the simulation clock. The new run starts
at zero, with empty history and the normal BDF1 startup before BDF2. Loading does
not perform a nonlinear solve or advance a timestep. Initialization belongs in
`initializeProblem()`, before assembly and time solver setup.

Format version, source physical time, fields, reference mesh/DOF fingerprints,
FE types, dimensions, components, global DOF counts and index bases are validated
before any solution block is replaced. Initial solution loading requires metadata;
`Allow legacy restart` does not bypass these checks. Old integration settings and
history are not imported, so the new dt, BDF order, solver and physical parameters
may differ. The caller supplies a suitable converged solution; loading does not
verify that it is an equilibrium for the new parameters.

`Initial solution` and `Restart` are mutually exclusive. FSI and Newmark initial
states are not supported yet because geometry, structural derivatives and outlet
state need additional initialization rules. Use a separate checkpoint output
directory when writing the new trajectory to preserve the source files.

[navierStokes_2D_3D_initialSolution](../feddlib/problems/tests/RestartTests/navierStokes_2D_3D_initialSolution/README.md)
generates converged steady BFS solutions in 2D and 3D on four ranks, then loads
them on four or six ranks. It checks exact initial field loading, zero time and
empty history, a zero-duration run, and the first two timesteps against direct
vector initialization, with different integration settings and no source history.

## Validation

With tests enabled and the required Trilinos solvers available:

```sh
cmake --build build --target problems_navierStokes_2D_3D_bfs problems_fsi_2D_turek problems_unsteadyNonLinElasticity_restart
ctest --test-dir build --output-on-failure -R 'problems_(navierStokes_2D_3D_bfs|fsi_2D_turek|unsteadyNonLinElasticity_restart)'
```

- `navierStokes_2D_3D_bfs` runs the 2D `BFS2d_1600.mesh` and 3D
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
- `fsi_2D_turek` runs an uninterrupted simulation and a restarted simulation on
  four MPI ranks, comparing fluid velocity, pressure, and solid displacement at
  the same final time with relative tolerance `1e-12`. Both phases load the case
  from `parametersProblemFSI.xml`; `parametersProblemFSI_restart.xml` contains
  only the overrides to resume at `0.01` and stop at `0.02`. The test disables
  visualization and benchmark exports and reports the relative error for each
  of the three fields.

- `unsteadyNonLinElasticity_restart` uses the 2D P2 `square_h02.mesh` unit square
  (3,015 vertices, 5,828 triangles), with the setup from the
  `unsteadyNonLinElasticity` example and Saint Venant-Kirchhoff
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
  Both rank pairs pass at `1e-12`; the largest relative error with `square_h02.mesh`
  was approximately `5.28e-14` (acceleration, `4 -> 6`). During initial validation
  on the coarse `square_solid.mesh`, scaling only velocity or only acceleration
  at `0.01` by `1.1` in separate disposable checkpoints made continuation fail
  the numerical comparison. The FSI regression also passed at that stage.

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
| 3: metadata and compatibility | Version-1 component manifests and shared pre-restoration validation, including global mesh/DOF fingerprints and HDF5 shape/type checks. | Positive restart and MPI rejection tests cover mesh/discretization/settings/history mismatches and explicit legacy mode. |
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

## Failure recovery

Nonlinear BDF, nonlinear Newmark and FSI runs can enable recovery independently
of ordinary checkpoint scheduling:

```xml
<ParameterList name="Timestepping Parameter">
  <Parameter name="Failure recovery" type="bool" value="true"/>
  <Parameter name="Recovery directory" type="string" value="recoveryCheckpoints"/>
  <Parameter name="Recovery linear iteration fraction" type="double" value="0.5"/>
  <Parameter name="Recovery interval steps" type="int" value="1"/>
</ParameterList>
```

The default fraction triggers when the mean linear iteration count of the
Newton solves in one timestep is strictly greater than half the active Belos
solver's `Maximum Iterations`. Direct solvers have no iteration-count trigger.
Zero Newton linear solves give a zero mean. Geometry solves are excluded.
Reaching `MaxNonLinIts` triggers only when the final permitted update has not
converged. Unsuccessful linear solves and nonfinite solution values also trigger
recovery.

`Cancel MaxNonLinIts`, in `Parameter`, controls stopping: when true, the driver
saves recovery information before raising a controlled failure; when false,
it warns and continues. A converged step with a high linear iteration count
continues with either setting. Ignoring nonconvergence permanently marks that
run's subsequent trajectory unreliable for recovery: later values cannot
replace the last reliable recovery state. Ordinary output retains its existing
behaviour. Recovery does not automatically retry or change the timestep.

The driver captures independent vector copies and full required history before
the solve. FSI captures coupled fields, fluid solution and mass-product history,
solid Newmark displacement/velocity/acceleration, ALE geometry and enabled outlet
history. Capture hooks retain the existing FSI/Newmark operation order and
start-of-step checkpoint layout. Consequently their latest completed capture
can precede the final solution by one timestep. Standalone BDF additionally
captures its final accepted state immediately.

Recovery generations are written to `generation_N.tmp`, with separate HDF5
handles, validated and closed before publication as `generation_N`.
`Complete.xml` marks a finished generation. `Latest.xml` selects the latest
complete generation and names the protected generation preceding the first
difficult timestep. It is replaced by a rename only after the new generation
is complete; older rolling generations are then removed. The protected
generation remains for the rest of the run. `RecoveryStatus.xml` records the
trigger, solver counts, convergence status, attempted time, cancellation decision
and preceding reliable time. Recovery files do not replace ordinary checkpoints.
Use a fresh recovery directory for each run; existing generations are never
truncated.

The recovery interval defaults to one step. Setting it to zero disables periodic
disk writes, retaining captures in memory until a trigger or an active warning
requires publication. A larger interval reduces routine disk output. A hard
process/node failure can only use generations already on disk; recovery does
not perform MPI/HDF5 operations in signal handlers.

To restart, read the directory and physical time from `Latest.xml` (or choose the
protected directory), then use the ordinary settings:

```xml
<Parameter name="Restart" type="bool" value="true"/>
<Parameter name="Restart directory" type="string" value="recoveryCheckpoints/generation_7"/>
<Parameter name="Time step" type="double" value="0.0175"/>
```

The usual metadata and compatibility validation applies. Restore normal solver
limits for continuation. Two four-rank 2D integration tests exercise warning,
stop and continue modes and compare recovery restarts against independently
computed reference solutions; a two-rank unit test checks criteria and interrupted
publication. Their setups are documented in
[fsi_2D_turek_recovery](../feddlib/problems/tests/RestartTests/fsi_2D_turek_recovery/README.md)
and [navierStokes_2D_3D_bfs_recovery](../feddlib/problems/tests/RestartTests/navierStokes_2D_3D_bfs_recovery/README.md).
