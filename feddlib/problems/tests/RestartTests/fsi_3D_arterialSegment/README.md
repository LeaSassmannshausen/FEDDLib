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
tolerance 1e-10 and flexible GMRES with tolerance 1e-12.

The `problems_fsi_3D_arterialSegment_Resistance`, `_Absorbing`, and `_Absorbing_Paper`
tests use separate working directories and four MPI ranks. Their transition time
is `0.00625`, between timesteps of size `0.0025`. Each test restarts from both `0.005`
and `0.01`, comparing velocity, pressure, and displacement at `0.015`, together
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

The variants also verify that changed averaging, missing outlet state and an
unknown state version are rejected. The Resistance variant additionally tests
`Safe all solution = true` with `Checkpointing = false`.

## Outlet checkpoint contract

An enabled pressure boundary model writes `FSIOutletState_<time>.xml` to the
configured checkpoint directory, alongside the field files. This versioned file
contains full-precision scalar history. The compatibility manifest identifies
the required state file and the pressure-model settings.

FSI retains its existing start-of-step checkpoint convention: a checkpoint at
`t_n` stores the solution at `t_n` and the outlet state **before** evaluating
the pressure load for the next solve. In particular, the retained previous flow
rate must not be replaced with a flow recomputed from the current velocity.
The uninterrupted run advances one extra step to write its final comparison
checkpoint; the restarted runs stop at the comparison time.

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
