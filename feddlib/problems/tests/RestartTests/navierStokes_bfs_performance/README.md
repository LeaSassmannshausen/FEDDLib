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
Results are in that directory's `results/timings.csv`. The `Performance` label
allows selection with `ctest -L Performance` or exclusion with `ctest -LE Performance`.

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
| `initial_solution_load` | Initial-field validation and loading, without time/history restoration |
| `paraview_resume` | Existing XMF/HDF5 validation, retained-frame selection and reopening |

The production hooks are enabled by `General/Checkpoint timings=true` and are
disabled by default elsewhere. They use Teuchos counters without adding MPI
barriers to the operations. `paraview_resume` covers exporter reopening; subsequent
frame writes and diagnostic text-log resumption are part of the phase's `advance`.

Standard output is CSV by default, also saved as `timings.csv`; `--verbose`
adds ordinary library progress output. Rows contain rank-minimum, mean and
**maximum wall time in seconds**, call counts, global/owned DOF counts, MPI size
and mesh settings. For scaling plots, use `seconds_max` as the parallel critical
path; divide by the reported call count when a per-call average is desired.
Each operation's row sums its calls within that phase. **Nested counters are
inclusive: do not add child times to their parent or to the phase total.**

Phase timers synchronize their starts; operation timers do not. `total` includes
setup, phase synchronization and exporter closure, but excludes CSV reporting.
Measurements use the actual filesystem with no cache flushing: repetitions and
later reads may benefit from filesystem caches. Inspect repetitions separately
and report warm/cold-cache conditions when publishing results.
