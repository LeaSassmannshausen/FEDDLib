# Checkpoint I/O optimization measurements

Measured locally on 2026-11-02. Each change was built and checked before the next
one was applied. Solver configuration and timestep settings were held fixed;
solver timings are outside this assessment.

## Measurement setup

Nine MPI ranks, 3D structured BFS, H/h=9, two producer steps, BDF2, three
repetitions and `Aspect=checkpoint_write`. The same executable/input interface
was used for the original baseline and all four cumulative stages:

```sh
mpirun -np 9 ./problems_navierStokes_bfs_performance.exe \
  --benchmarkfile=parametersBenchmark.xml \
  --dimension=3 --H-h=9 --steps=2 --restart-steps=1 --repetitions=3 \
  --aspect=checkpoint_write --output-directory=NEW_OUTPUT_DIRECTORY
```

Every run used a fresh output directory on the local `/tmp` filesystem. These
measurements do not establish performance on the cluster's Lustre filesystem.
No cache flushing or minimum-speedup acceptance threshold was used.

## Cumulative changes

1. **Replicated vector properties.** In `HDF5Toolbox`, integer, real and string
   datasets are created collectively and written collectively with only rank
   zero selecting data. Other ranks select no elements. The existing names,
   types and shapes of `GlobalLength`, `NumVectors` and `__type__` are preserved.
   Preparation/write failures are reported collectively.
2. **Internal HDF5 metadata.** `HDF5Toolbox::create()` enables
   `H5Pset_coll_metadata_write` on the MPI file-access property list. This controls
   HDF5's internal metadata cache, independently of the application properties in
   stage 1. Collective metadata reads were not enabled.
3. **Batch flushing.** `TimeProblem` writes current/history fields with
   `flushAfterWrite=false`, then checks one flush per field file before
   publishing checkpoint metadata. Its accepted-step and legacy FSI/Newmark
   export paths both use this completion boundary. Recovery snapshots batch each
   file and check its flush on explicit closure before persisted-data validation
   and `Complete.xml` publication. Standalone exporter calls still flush by
   default. The former second flush after every standalone write was removed.
   Failed vector-dataset creation is rejected collectively before property writes.
4. **Write-buffer/importer reuse.** `HDF5Toolbox` caches its contiguous vector
   buffer and importer using the vector's actual Tpetra storage map. It rebuilds
   when that map's identity or the column count changes and releases the cache
   on closure. Rebuild decisions are collective; partial-rank replacement is
   rejected before any dataset is created. Importer validation in
   `MultiVector::importFromVector()` uses actual storage maps, matching importer
   construction, rather than potentially copied FEDD map wrappers.

## Local timing results

All times below are seconds. The mean/median/range are across the three reported
MPI mean values. The MPI-max column averages the three reported MPI maximum
values. Each `checkpoint_write` includes both BDF2 states, velocity/pressure,
the required flushes and checkpoint metadata. File closure is measured separately.

| Stage | Checkpoint mean | Checkpoint median | Repetition range | Mean MPI maximum |
|---|---:|---:|---:|---:|
| Original baseline | 0.139981 | 0.096031 | 0.085564–0.238348 | 0.140004 |
| 1: rank-zero property payloads | 0.163425 | 0.095223 | 0.093170–0.301882 | 0.163576 |
| 2: collective internal metadata | 0.110143 | 0.110962 | 0.095091–0.124377 | 0.110351 |
| 3: batched checked flushing | 0.071244 | 0.070176 | 0.063492–0.080066 | 0.071263 |
| 4: buffer/importer reuse | 0.087959 | 0.051249 | 0.050172–0.162455 | 0.087967 |

Component means below are milliseconds. Flush time sums the immediate and
explicit flush counters where both existed; the counters measure disjoint calls.

| Stage | Vector properties | Redistribution | HDF5 flushes | Explicit flush calls/checkpoint |
|---|---:|---:|---:|---:|
| Original baseline | 2.242 | 29.400 | 49.854 | 8 |
| 1 | 8.073 | 47.682 | 51.041 | 8 |
| 2 | 2.602 | 38.092 | 51.869 | 8 |
| 3 | 2.574 | 27.051 | 26.374 | 2 |
| 4 | 2.741 | 15.102 | 22.903 | 2 |

The nine-rank run does not demonstrate a property-write speedup from stage 1;
its collective coordination has a cost. Stage 2's isolated benefit is also
unresolved by these local samples. Both changes need cluster measurements.

Stage 3 verifies the intended reduction from eight explicit flush calls to two.
Stage 4 verifies two buffer/importer setups for four field writes, reusing each
field's buffer for its history state. Its redistribution timer falls from
27.051 ms after stage 3 to 15.102 ms. One stage-4 repetition spent 95.645 ms
creating files, so its checkpoint mean increased despite lower redistribution
and flushing costs. The ranges and medians retain that variation.

`fields_flush` is now a sibling of `fields_write` beneath `export` and contains
`hdf5_flush`. `hdf5_flush_after_write` is absent from batched Navier-Stokes writes.
Compare `checkpoint_write` totals across stages; `fields_write` alone no longer
includes flushing. `hdf5_buffer_setup` is a child of `hdf5_redistribution` and
includes initial buffer/importer preparation and the first import.

## Correctness checks

After every stage, all six benchmark HDF5 files were compared with the baseline:
dataset keys, dimensions, types and values were identical, with zero numerical
difference. The file format and restart-reader interface were preserved.

The new `problems_tests_unitTests_hdf5Checkpoint_MPI_2` checks rank-zero scalar
payloads, empty writer/reader ranks, changed read distribution, deferred and
immediate flushing, explicit closure with pending data, closed-file rejection,
copied-vector cache reuse, collective map replacement, ownership/local ordering,
column count, global length and index base. It also verifies collective rejection
of partial-rank map replacement before I/O.

The final integration selection covers metadata, recovery and output-history
unit tests plus all restart tests, including processor changes, FSI boundary
conditions, Newmark, recovery, initial solutions, variable timesteps and the
three performance smoke tests:

**Result: all 27 tests passed in 594.25 seconds.** The final performance reports
also verify the timer hierarchy and two buffer setups/two explicit flushes for
four field writes in every measured checkpoint branch.

```sh
ctest --test-dir build \
  -R '^problems_tests_(restartTests_|unitTests_(hdf5Checkpoint|checkpointMetadata|recoveryCheckpoint|outputHistory)_MPI_2$)' \
  --output-on-failure
```

Raw reports, checkpoint files, JSON summaries, build logs and cumulative source
patches for each stage are retained at:
`/tmp/feddlib-checkpoint-optimizations-20261102-chn9cge9`.

For the cluster comparison, retain the original input files, aspect, H/h,
timestep counts, placement and filesystem settings at 243, 576 and 1125 ranks.
Rebuild the executable and use a fresh output directory for every run. Inspect
`checkpoint_write`, `hdf5_vector_metadata_write`, `fields_flush`,
`hdf5_redistribution` and the buffer-setup counts; report each repetition as well
as aggregate statistics.
