# FSI restart tests

The driver contains only the 3D arterial segment: P1/P1 fluid, a P1 Neo-Hooke
wall, explicit ALE geometry, and the FaCSI preconditioner. The supplied meshes
are used without scaling. In centimetre units their length is 0.05 cm (0.5 mm)
and their fluid radius is 0.09 cm.

A parabolic inlet profile is computed from mesh coordinates and normalized by
BCBuilder to the prescribed volumetric flow. The old stored profile files belong
to a different mesh and are not used. The supplied material properties and flow
amplitude are retained; ramp/transition times and the simulation duration are
shortened for this regression test. The solver uses Newton with relative residual
tolerance 1e-10 and flexible GMRES with tolerance 1e-10.

The `problems_fsi_3D_arterialSegment_Resistance`, `_Absorbing`, and `_Absorbing_Paper`
tests use separate working directories and four MPI ranks. Each loads frozen
references and advances **exactly two timesteps per simulation**, with
`dt = 0.0025`. Resistance runs once, from `0.005` to `0.01`, while writing
save-all checkpoints. Each absorbing variant runs from both `0.005` to `0.01`
and `0.01` to `0.015`, before and after its `0.00625` transition time.
The tests compare velocity, pressure and displacement at the corresponding endpoint, together
with initial areas, transition area, flow-rate history, outlet pressure and flags.
The field comparison uses the combined l2 bound
`absolute error <= 1e-15 + 1e-12 * reference l2 norm`. The absolute tolerance
is in the stored field's units and covers round-off for small displacements;
both bounds are configurable in the timestepping parameters. Scalar history uses
`absolute error <= 1e-14 + 1e-12 * abs(reference)`. Solver output is streamed to
CTest and can be viewed with `-V`.

The original 2D Turek restart test remains in the sibling `fsi_2D_turek` folder.

The resistance variant first checks the integrated traction for an affine
velocity against its analytic value, covering pressure and viscous components.
This check uses a non-unit density: the input viscosity is kinematic, and the
viscous traction uses dynamic viscosity (density times kinematic viscosity).

Normal CTest runs perform no uninterrupted reference simulation. Complete
checkpoint bundles are stored in
`UnitTest/ReferenceSolutions/fsi_3D_arterialSegment/{Resistance,Absorbing,Absorbing_Paper}`.
The build copies them into the test's binary directory. Every CTest case makes
its own writable copy; mutation/rejection checks leave the frozen bundles intact.
The driver verifies the completed timestep-counter difference, in addition to
comparing fields and scalar history. The runner explicitly passes
`--restart-steps=2`. For manual runs, a positive value overrides the stopping
time with a fixed step count; zero (the default) preserves the XML `Final time`.
This keeps the separate variable-timestep regression's configured endpoints.

The variants also verify that changed averaging, missing outlet state and an
unknown state version are rejected. These are validation-only launches without
timestepping. The single Resistance simulation uses `Safe all solution = true`
with `Checkpointing = false` and writes every solution in `saveAllCheckpoints`.
The runner requires outlet-state files at `0.005` and `0.0075`, without another
simulation. Its analytic traction check and field/scalar error bounds are retained. Each
simulation prints a phase description before the solver output.

## Outlet checkpoint contract

An enabled pressure boundary model writes `FSIOutletState_<time>.xml` to the
configured checkpoint directory, alongside the field files. This versioned file
contains full-precision scalar history. The compatibility manifest identifies
the required state file and the pressure-model settings.

FSI retains its existing start-of-step checkpoint convention: a checkpoint at
`t_n` stores the solution at `t_n` and the outlet state **before** evaluating
the pressure load for the next solve. In particular, the retained previous flow
rate must not be replaced with a flow recomputed from the current velocity.
Only explicit reference generation advances one extra step to write its final
comparison checkpoint; normal tests stop at their two-step comparison endpoints.

Restoration validates the outlet snapshot before importing any coupled solution
or ALE field. Rank zero reads the scalar state and broadcasts it to all ranks.
Initial areas and captured transition areas are retained; the pressure model is
not advanced during restoration. A pending transition is captured once on the
first pressure evaluation reaching its time.

Old field-only checkpoints remain usable for the `None` model under the existing
compatibility rules. An enabled pressure model requires its outlet snapshot,
including when `Allow legacy restart` is enabled.

Run these tests with:

```sh
ctest --test-dir build -V -R 'problems_fsi_3D_arterialSegment'
```

The metadata unit test additionally checks scalar round-trip precision, invalid
areas, missing history, nonfinite values, missing transition state, and rejection
of changed density or kinematic viscosity for the resistance model.

## Deliberate reference regeneration

After reviewing changes to the geometry, physics, integration or solver settings,
build the driver and run the generated maintenance script into a new directory:

```sh
cmake \
  -DEXECUTABLE="$PWD/build/feddlib/problems/tests/RestartTests/fsi_3D_arterialSegment/problems_fsi_3D_arterialSegment.exe" \
  -DREFERENCE_OUTPUT_DIRECTORY=/tmp/new-arterial-references \
  -P build/feddlib/problems/tests/RestartTests/fsi_3D_arterialSegment/generateOutletReferences.cmake
```

This runs each model uninterrupted on four ranks to `0.0175`, with `dt=0.0025`.
Resistance uses save-all output; the absorbing models save checkpoints at
`0.005`, `0.01` and `0.015` with their required prior history. The generated model
directories contain field/ALE/derivative/RHS HDF5 files, manifests, scalar outlet
snapshots and generation notes. Review them before replacing the corresponding
source reference directories; rebuild to refresh the binary copies. The script
rejects an existing output destination and never replaces committed references
automatically. Normal CTest runs never invoke it.
