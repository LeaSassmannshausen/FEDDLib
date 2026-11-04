# fsi_artery_boundary_tests

Independent boundary-condition benchmark copied from `fsi_artery`. The driver,
solver settings, preconditioner settings, and problem XML can be edited here.
The initial configuration reuses the existing artery meshes and inlet profiles.

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
2. Supply inlet profiles matching the geometry, discretization, and node
   numbering. The driver reads dataset `solution` from
   `<Inflow Profile Prefix><Discretization>.h5`; set `Parameter Fluid/Inflow
   Profile Prefix` accordingly. Add the profiles to the CMake copy list or use
   a path prefix. The inherited profile must be normalized to one.
3. Match the boundary labels in `main.cpp` to the new meshes. The inherited
   labels are fluid inlet **4**, fluid outlet **5**, coupled interface **6**
   with inlet/outlet rings **9/10**, solid inlet/outlet **7/8**, and solid
   support strips **13/14**. Review the fluid, structure, and ALE boundary
   sections together when changing them.
4. Review the inherited `readAndPartition(15, "mm", true)` call: the original
   case converts millimetre mesh coordinates to centimetres. Material, fluid,
   and flow-rate parameters must use consistent units. The inherited inflow
   callback and structural constraints also assume their original directions.

The copied problem defaults to a 3D P1 Neo-Hookean solid. For AceGen SCI
materials, use P2 and select the intended material in `Parameter Solid`;
see [the AceGen FSI notes](../../../../../doc/acegen-fsi.md).
