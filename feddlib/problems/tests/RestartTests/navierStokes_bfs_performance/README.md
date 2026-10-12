# Structured BFS checkpoint performance

This executable measures the production Navier-Stokes checkpoint/restart paths
on structured P1-P1 backward-facing-step meshes, in 2D or 3D with BDF2. There
are no reference-solution comparisons or runtime pass/fail thresholds. Solver
nonconvergence, invalid input and I/O failures still fail the run.

## Local tests

Build `problems_navierStokes_bfs_performance`, then run:

```sh
ctest --test-dir build -R '^problems_navierStokes_bfs_performance_' -V
```

Two short cases use nine MPI processes, H/h=3, four producer steps and two
continuation steps. CTest recreates only each case's dedicated working directory.
Results are in that directory's `results/timings.csv` and `results/stacked_timings.txt`. The `Performance` label
allows selection with `ctest -L Performance` or exclusion with `ctest -LE Performance`.
A third nine-process case verifies the parameter-file interface: two 2D restart
samples with H/h=5, launching from outside the input-file directory.

## Configure a run with a parameter file

Edit [`parametersBenchmark.xml`](parametersBenchmark.xml), then launch the executable:

```sh
mpirun -np 9 ./problems_navierStokes_bfs_performance.exe \
  --benchmarkfile=parametersBenchmark.xml
```

The supplied file selects `Aspect=restart`, dimension 2, H/h=5, two producer
steps, one continuation step and two independent repetitions. Its output goes
to `performance_restart` in the directory from which you launch the executable.
Each repetition first creates its own checkpoint, then measures restoration and
continuation. This runs directly on nine processes without editing other files.

The controls are distinct from `parametersProblem.xml`, which sets physics, dt
and nonlinear convergence settings. The benchmark XML references the problem,
preconditioner and linear-solver files; those paths resolve relative to the
benchmark XML's directory. Thus an absolute benchmark-file path also works when
launching from a different shared-scratch working directory. Output paths resolve
relative to the launch directory, or can be absolute shared-filesystem paths.

| Benchmark parameter | Meaning |
|---|---|
| `Dimension` | 2D or 3D structured BFS |
| `H/h` | Mesh subdivisions per subdomain edge, at least 2; fixes the local element workload |
| `BFS length` | Downstream length; also determines supported rank counts |
| `Aspect` | Select the measurement scenario from the table below |
| `Steps` | Producer/baseline step count; final checkpoint time is `Steps × dt` |
| `Restart steps` | Continuation length, or new steps from time zero for initial fields |
| `Repetitions` | Independent samples with separate files and timer branches |
| `Output directory` | New destination for timing reports, effective inputs and checkpoints |
| `Quiet`, `Stacked timers` | Control solver progress and the additional hierarchical report |

| `Aspect` value | Phases executed | Timings to inspect |
|---|---|---|
| `baseline` | `baseline` only | `spatial_assembly`, `advance`, `total` without checkpoint/export I/O |
| `checkpoint_write` | `checkpoint_run` only, without visualization/text exports | `checkpoint_write`, `fields_write`, `metadata_write` |
| `restart` | `checkpoint_run`, `restart` | `compatibility_validate`, `primary_fields_read`, `bdf_history_read` |
| `initial_solution` | `checkpoint_run`, `initial_solution` | `initial_solution_load`, `primary_fields_read`; no history restoration |
| `resume_output` | `checkpoint_run` with exports, `resume_output` | `paraview_resume`, and `advance` for resumed trajectory/export costs |
| `all` | Existing complete benchmark suite | All stages and operations |

Prerequisite producer time is reported separately from the selected aspect.
For example, when assessing restart costs, read the `restart` rows in the CSV
or the `restart` branch in the stacked report. No solutions are compared.
You cannot infer isolated restart costs from the complete invocation time,
because that also includes the producer, solves, setup and timing-report I/O.

For a 3D cluster series, for example, set `Dimension=3`, `H/h=8`, `BFS length=4`,
`Aspect=restart`, `Steps=4` and `Repetitions=3`, then launch on 9, 72, 243, ...
ranks. Keep these settings and the physics/solver files fixed across the series;
choose a different output directory for every allocation. MPI process count is
set by the launcher, not by the XML. This example uses the default one-level
preconditioner; select an appropriate cluster preconditioner in the referenced
file when needed, and use the same one throughout a series.

Explicit command-line options override XML controls. For example, to run that
same file on 2D with a different output location:

```sh
mpirun -np 36 ./problems_navierStokes_bfs_performance.exe \
  --benchmarkfile=parametersBenchmark.xml \
  --dimension=2 --H-h=12 --output-directory=weak_2D_36
```

Invocations without `--benchmarkfile` load `parametersBenchmark.xml` by default.
Select `--aspect=all` to run the full suite; the local full-suite CTest cases do
this explicitly and override the mesh/step/repetition controls for short runs.
`--aspect=restart` also selects a single aspect without a benchmark XML. For a
single selected aspect, its phase selection determines which scenarios run;
the older `--no-baseline`, `--no-resume-output`, `--no-initial-solution` switches
customize `Aspect=all` only. `parametersEffective.xml` records the resolved
controls, MPI layout and input-file paths; phase XML files retain that record
alongside each phase's checkpoint/time settings.

## Larger runs and weak scaling

From the built executable/input directory, for example:

```sh
mpirun -np 36 ./problems_navierStokes_bfs_performance.exe \
  --dimension=2 --H-h=12 --repetitions=3 --output-directory=weak_2D_36

mpirun -np 72 ./problems_navierStokes_bfs_performance.exe \
  --dimension=3 --H-h=8 --repetitions=3 --output-directory=weak_3D_72
```

Use your cluster's MPI launcher/allocation. No rank count is compiled into the
executable: it reads the MPI communicator size. `--problemfile`, `--precfile`
and `--solverfile` accept alternate XML files. `--help` lists all options.
Choose a new output directory for every invocation: existing directories are
rejected to preserve prior results. Each repetition has independent checkpoint
files and preserves its effective phase parameter files.

For the default downstream length of 4, the existing structured BFS generator
requires **ranks = 9 × k^dimension**, with an integer k >= 1:

| Dimension | Supported rank counts |
|---|---|
| 2D | 9, 36, 81, 144, 225, 324, ... |
| 3D | 9, 72, 243, 576, 1125, 1944, ... |

For other `--length=L`, ranks must be `(2L+1) × k^dimension`. Unsupported layouts
are rejected before building a mesh. Keep dimension, length, **H/h**, timestep
count, dt, solver settings and export cadence fixed across a weak-scaling series.
The driver sets N=k and the mesh spacing h=1/(k × H/h), retaining the same number
of elements per rank (2 × (H/h)^2 triangles or 6 × (H/h)^3 tetrahedra). The physical
BFS stays fixed as it is refined. This is computational weak scaling; changing
h can still affect solver iterations, which is why solve and I/O times are separate.
Owned DOFs differ slightly at boundaries; CSV includes their rank minimum/maximum.

The default preconditioner is a one-level FROSch configuration for short local
runs. For large jobs, select an appropriate scalable preconditioner with
`--precfile`; keep it identical across the scaling series. Use a consistent
optimized build, MPI placement, thread count and filesystem for measurements.

## Phases and timing rows

The default invocation executes:

1. `baseline`: fresh solve with checkpoint and visualization/text output disabled.
2. `checkpoint_run`: fresh solve, writes the final accepted checkpoint with BDF2
   history and produces visualization/text output for the continuation phase.
3. `restart`: restores that checkpoint and advances with checkpoint/output I/O disabled.
4. `resume_output`: independently restores the same checkpoint, continues the
   producer's ParaView/text files and writes a new final checkpoint in a separate directory.
5. `initial_solution`: loads the checkpoint's primary fields at time zero and
   advances a fresh trajectory without restoring the old time or BDF history.

`--steps` controls producer/baseline length; `--restart-steps` controls the other
phases. `--no-baseline`, `--no-resume-output` and `--no-initial-solution` skip their
respective optional phases. With `--no-resume-output`, the producer also disables
visualization/text output, providing an isolated checkpoint-write measurement.
Recovery/failure injection and FSI/Newmark are outside this Navier-Stokes benchmark.

Each phase reports `mesh_build`, `problem_initialization`, `spatial_assembly`,
`time_setup`, `advance` and `total`, plus the operations actually executed:

| Operation | Timed work |
|---|---|
| `metadata_prepare` | Mesh/DOF fingerprints and schema preparation; includes validation on restart |
| `compatibility_validate` | Manifest, integration/field/mesh and required-history compatibility checks |
| `primary_fields_read` | Collective HDF5 velocity/pressure import and installation |
| `bdf_history_read` | Collective HDF5 reads and allocation of BDF solution history |
| `fields_write` | Collective HDF5 field writes, including exporter creation on the first write |
| `metadata_write` | Checkpoint manifest serialization/publication |
| `checkpoint_write` | Complete accepted-state checkpoint, including history and metadata |
| `check_for_export` | Check save-all/checkpoint flags, initialize the schedule if needed and decide whether output is due; also measured when no write occurs |
| `export` | Actual checkpoint export of fields/history and metadata, for accepted-step Navier-Stokes and the FSI/Newmark path |
| `check_for_export_and_export` | Entire FSI/Newmark `TimeProblem::checkForExportAndExport` call, including optional recovery capture, decision checks and any export |
| `initial_solution_load` | Initial-field validation and loading, without time/history restoration |
| `paraview_resume` | Existing XMF/HDF5 validation, retained-frame selection and reopening |

The production hooks are enabled by `General/Checkpoint timings=true` and are
disabled by default elsewhere. They use Teuchos counters without adding MPI
barriers to the operations. `paraview_resume` covers exporter reopening; subsequent
frame writes and diagnostic text-log resumption are part of the phase's `advance`.
`export` measures checkpoint output, rather than ParaView frame output. In the
Navier-Stokes accepted-step path it is nested inside `checkpoint_write`; the
`check_for_export` decision is timed separately before that write. FSI/Newmark
uses `check_for_export_and_export` as the parent, with separate decision and
write regions. Its `export` region also contains `fields_write` and
`metadata_write`. The full legacy parent includes recovery snapshot capture
when enabled; that capture is outside the decision/write regions. Its timers
appear in stacked reports when those paths are exercised, rather than in this
standalone Navier-Stokes run. `fields_write` counts one call per accepted-step
state, or per legacy field/history batch, so compare call counts accordingly.

Standard output contains CSV timing rows followed by a stacked timer report.
The CSV is also saved separately as `timings.csv`; `--verbose` adds ordinary
library progress output. Rows contain rank-minimum, mean and
**maximum wall time in seconds**, call counts, global/owned DOF counts, MPI size
and mesh settings. For scaling plots, use `seconds_max` as the parallel critical
path; divide by the reported call count when a per-call average is desired.
Each operation's row sums its calls within that phase. **Nested counters are
inclusive: do not add child times to their parent or to the phase total.**

Stacked timing is enabled by default (`--stacked-timers`) and saved as
`stacked_timings.txt` beside the CSV. `--no-stacked-timers` disables the stack
and keeps CSV-only output. The hierarchy separates each repetition and phase:

```text
BFS restart performance
  repeat_1
    checkpoint_run
      advance
        FEDD checkpoint - checkpoint_write
          FEDD checkpoint - export
            FEDD checkpoint - fields_write
            FEDD checkpoint - metadata_write
    restart
      problem_initialization
        FEDD checkpoint - metadata_prepare
          FEDD checkpoint - compatibility_validate
        FEDD checkpoint - primary_fields_read
      advance
        FEDD checkpoint - bdf_history_read
```

Existing solver/preconditioner `TimeMonitor` counters appear automatically
under the stages that execute them. This requires Trilinos's
`Teuchos_ENABLE_STACKED_TIMER_IN_TIME_MONITOR=ON` (the default).
The stacked report provides MPI mean,
minimum/maximum times, the ranks attaining those extrema, call counts and
fractions of parent time. `Remainder` denotes time outside the measured children;
the root/repetition regions include bookkeeping and CSV reporting (`csv_report`).
These parent fractions describe inclusive timings, not extra costs to add.
Stages use the same synchronization policy as the CSV; attaching the stacked
timer adds no barriers to production checkpoint operations. A detected library
timer overlap fails the run rather than publishing an incomplete stack.

Phase timers synchronize their starts; operation timers do not. `total` includes
setup, phase synchronization and exporter closure, but excludes CSV reporting.
Measurements use the actual filesystem with no cache flushing: repetitions and
later reads may benefit from filesystem caches. Inspect repetitions separately
and report warm/cold-cache conditions when publishing results.
