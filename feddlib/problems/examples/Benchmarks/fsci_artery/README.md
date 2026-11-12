# fsci_artery

Fluid–solid–chemistry interaction benchmark. The supplied FSCI problem replaces
FSI's solid equation with the current SCI solid/chemistry system. Chemistry is
solved in the artery wall; this case does not add a fluid transport equation.

## Build and run

From the workspace root:

```sh
cmake -S source/FEDDLib -B build
cmake --build build --target problems_fsci_artery -j 4
cd build/feddlib/problems/examples/Benchmarks/fsci_artery
mpirun -np 4 ./problems_fsci_artery.exe --problemfile=parametersProblemFSCI.xml
```

This is a manually run benchmark, like `fsi_artery_boundary_tests`; there is no
frozen CTest reference solution. The default is the same 12 mm fluid/solid pair,
3D P2 displacement, velocity and concentration, and P1 fluid pressure. CMake
copies all mesh pairs and P2 inlet profiles supported by that FSI benchmark,
including the original 8 mm pair. Choose both `Mesh Partitioner/Mesh 1 Name`
(fluid) and `Mesh 2 Name` (solid) as described in its README. The driver selects
`laplace_parabolic_<fluid mesh basename>_P2.h5` automatically. Profiles remain
shared source assets in `fsi_artery_boundary_tests`, rather than duplicate HDF5
files here.

The fluid inlet is flag 4 at z = 0, directed along +z; outlet flag 5 and the
existing ramped `Parameter Fluid/Flowrate` normalization are retained. The solid
supports and ALE boundary conditions match the FSI benchmark.

## Chemistry and preconditioning

`Parameter/FSCI = true` selects the coupled time integration path. Geometry may
be implicit or explicit through `Parameter/Geometry Explicit`. With
`Parameter/Chemistry Explicit = false`, concentration joins the monolithic
system after the interface (GE) or geometry (GI) block. With chemistry explicit,
SCI advances transport separately using its existing staggered scheme and a
separate transport mesh.

`General/Preconditioner Method = FaCSCI` uses the supplied extension of FaCSI:
the coupled solid/chemistry inverse comes from `parametersPrecSCI.xml`, fluid
settings from `parametersPrecFluidMono.xml`, and geometry settings from
`parametersPrecGeometry.xml`. Staggered transport uses `parametersPrecChem.xml`.
The `precfileSCI` and `precfileChem` command-line arguments override these files.

`Parameter Diffusion` controls the concentration boundary: `Inflow Flag = 11`
selects the outer wall, `Inflow Start Time = 0` and `Inflow Concentration = 1`
prescribe the default concentration. The selected AceGen SCI material provides
the coupled chemical residual and tangent; its material diffusivity is set in
`Parameter Solid/1/D0`. Use consistent units and material parameters, as in the
standalone SCI case. `Parameter Diffusion/D0` supplies the separate transport
operator for staggered chemistry.

Fluid uses the existing BDF integration, solid uses Newmark, and implicit SCI
chemistry retains its native backward-Euler concentration rate. Material history
advances once after each accepted coupled step. Existing `Exporter` flags for
material fields, element stresses, Gauss-point stresses and P1 stress recovery
also apply to FSCI.

SCI's existing restrictions carry over: the removed linear SCI option and legacy
material-history file I/O/restart are unavailable. This benchmark must start
from an initial state with checkpointing/restart disabled.

## Validation and open issue

The executable builds. Small-mesh, single-step checks cover implicit and explicit
geometry, coupled and staggered chemistry, and a two-rank NOX run. Nodal,
element-average, Gauss-point and P1 stress outputs were checked for finite data.
The existing four-rank FSI regression also passes.

Longer small-mesh runs still lose linear/nonlinear convergence after the first
update or timestep. A comparison with the existing FSI driver also showed
convergence trouble beyond its first step; the cause has not been isolated.
These checks do not establish convergence of the default artery simulation.
The benchmark stops on nonlinear failure rather than exporting an unconverged
state. Material algorithms and the existing FSI solver were not changed to
address this separate numerical issue.
