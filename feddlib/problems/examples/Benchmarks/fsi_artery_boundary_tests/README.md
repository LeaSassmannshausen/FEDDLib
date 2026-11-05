# fsi_artery_boundary_tests

Independent boundary-condition benchmark copied from `fsi_artery`. The driver,
solver settings, preconditioner settings, and problem XML can be edited here.
The source defaults reuse the existing artery meshes and prescribe an analytic
elliptical inlet profile. No HDF5 inlet profile is loaded.

## Build and run

From the workspace root:

```sh
cmake -S source/FEDDLib -B build
cmake --build build --target problems_fsi_artery_boundary_tests -j 4
cd build/feddlib/problems/examples/Benchmarks/fsi_artery_boundary_tests
mpirun -np 4 ./problems_fsi_artery_boundary_tests.exe --problemfile=parametersProblemFSI.xml
```

The configure command assumes an existing configured FEDDLib build. This is a
manually run benchmark; it has no CTest registration or reference solution yet.

## Add the geometry

1. Add the fluid and solid `.mesh` files and update the mesh-copy list in
   `CMakeLists.txt`. Set `Mesh Partitioner/Mesh 1 Name` and `Mesh 2 Name` in
   `parametersProblemFSI.xml` to the new filenames. Absolute paths in the XML
   can also point to meshes outside the build directory.
2. Set the analytic inlet parameters in `Parameter Fluid`: `BC Type` must be
   `parabolic_artery`. Configure `Flowrate` (cm³/s), `Inflow Ramp Time`, `Inflow
   Center X/Y/Z`, `Inflow Direction X/Y/Z`, `Inflow Axis X/Y/Z`, and `Inflow
   Radius 1/2`. These describe a centre, inward normal, in-plane principal
   axis, and two ellipse radii. Lengths are in **centimetres after conversion**.
3. Match the boundary labels in `main.cpp` to the new meshes. The inherited
   labels are fluid inlet **4**, fluid outlet **5**, coupled interface **6**
   with inlet/outlet rings **9/10**, solid inlet/outlet **7/8**, and solid
   support strips **13/14**. Review the fluid, structure, and ALE boundary
   sections together when changing them.
4. Review the inherited `readAndPartition(15, "mm", true)` call: the original
   case converts millimetre mesh coordinates to centimetres. Material, fluid,
   and flow-rate parameters must use consistent units. The inherited inflow
   structural constraints also assume their original directions.

The copied problem defaults to a 3D P1 Neo-Hookean solid. For AceGen SCI
materials, use P2 and select the intended material in `Parameter Solid`;
see [the AceGen FSI notes](../../../../../doc/acegen-fsi.md).

## Analytic inlet

Only `parabolic_artery` is supported by this benchmark. It sets
`u = amplitude(t) * max(0, 1 - (s/a)^2 - (q/b)^2) * direction`,
where `s` and `q` are coordinates in the inlet plane. Direction and axis vectors
are normalized; the axis is first projected onto the inlet plane. The existing
flow-rate boundary-condition mechanism integrates the unit vector profile on
the current inlet mesh and sets `amplitude = Flowrate * ramp(t) / abs(profile flux)`.
The integration uses finite-element nodal values, leaves non-inlet nodes
(including coupled rings) zero in this profile, and sums contributions across MPI ranks. It is
repeated when applying the boundary values, so mesh deformation is accounted for.
A positive `Inflow Ramp Time` applies a half-cosine ramp to the requested flow
rate; zero disables the ramp. No HDF5 inlet profile is required. `MaxVelocity`
is no longer used. `Flowrate` also remains an input to the absorbing outlet model.

This normalizes the prescribed profile's contribution to the inlet flux.
The rings retain their FSI coupling; any normal velocity there adds to the
total flux. Enforcing the target including that unknown contribution would
require an additional coupling change.

The same analytic Dirichlet condition is registered for the coupled problem
and its fluid component. The fluid wall is the coupled solid interface, so no
zero wall-velocity condition is added there. `Inflow Flag` defaults to 4 and
must match `General/Flag Inlet Fluid`; the outlet and structural conditions
keep their existing configuration.

The active build-directory case is also saved as
`parametersProblemFSI_parabolic_artery.xml`. It uses the new fluid and solid
meshes already present there. Its inlet at `z=0` points along `+z`; the fitted
centre is approximately `(0.00363215, 0.00843393, 0)` cm and the radii are
`0.14054931` and `0.12337601` cm. The source default ellipse is fitted to the
original plaque geometry instead. Update these entries when changing meshes.

To retain build-directory mesh settings across later builds, run with
`--problemfile=parametersProblemFSI_parabolic_artery.xml`.
