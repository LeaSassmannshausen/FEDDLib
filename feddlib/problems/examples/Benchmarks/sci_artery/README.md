# sci_artery

Solid–chemistry (pharmaco-mechanical) artery benchmark. It uses the solid mesh,
material parameters and displacement constraints from `fsi_artery_boundary_tests`.
The coupled unknowns are displacement `d_s` (block 0) and concentration `c`
(block 1). It contains no fluid, FSI interface multiplier or ALE subproblem.

## Build and run

From the workspace root, with an existing configured build and AceGen Interface:

```sh
cmake -S source/FEDDLib -B build
cmake --build build --target problems_sci_artery -j 4
cd build/feddlib/problems/examples/Benchmarks/sci_artery
mpirun -np 4 ./problems_sci_artery.exe --problemfile=parametersProblemSCI.xml
```

This is a manually run benchmark. CMake copies its XML settings and the 10 mm
and 12 mm solid meshes from `meshes/artery_adaptations`. The default is 12 mm;
change `Mesh Partitioner/Mesh 1 Name` to the `_10mm.mesh` variant to use 10 mm.
Mesh coordinates are converted from mm to cm
using the same `readAndPartition(15, "mm", true)` call as the FSI case.

## Configuration

The default is a 3D P2/P2 nonlinear SCI problem using
`SCI_SMC_Active_Growth_Reorientation`, with the FSI benchmark's `Parameter Solid`
material sublist copied unchanged. Active response, growth and reorientation
switches therefore retain that case's settings. `SCI_NH` and `SCI_SMC_MLCK`
can be selected with matching material settings; the linear SCI option is removed.

The default mechanical load replaces fluid traction with a prescribed pressure
on the inner wall (flag 6). `Parameter/Volume force` gives its magnitude and
`Load Ramp End` the ramp duration, following the supplied SCI example's pressure
callback. Solid constraints use flags 7–14 exactly as in the FSI benchmark.
Chemistry prescribes `Inflow Concentration` on `Chemistry Flag` (default 6), from
`Inflow Start Time` onward. Other chemical boundaries have natural conditions.
Lengths, material data, pressure and diffusivity use the cm–kg–s convention.

`Parameter/SCI = true` selects the dedicated `advanceInTimeSCI` path.
The supplied time-loop ordering is retained, with Newmark mechanics and the
material's backward-Euler concentration rate. `Timestepping Parameter/BDF` is 1.
`Load Step Size` must match the initial `dt`, as required by the supplied loop.
`Load Stepping = true` retains its quasi-static option with zero solid inertia.
`Chemistry Explicit` retains the supplied staggered chemistry option; the default
is the coupled implicit system.

`General/ParaViewExport` exports displacement and concentration. Set
`Exporter/Export material fields` (or legacy `General/Export Stress`) to true for
the current 17-field material layout, including stress, energy, fibre directions,
stretches and active-state variables. Material fields are written to
`solid_material.h5` / `.xmf`.

Legacy material-history file import/export and SCI restart/checkpointing are
unavailable in this port and produce a clear error if enabled. Material history
still advances in memory at each SCI step. Other solvers' restart paths are
unchanged. The source files and time-loop body are based on the supplied
Hackathon SCI implementation, with API compatibility changes for this branch.
