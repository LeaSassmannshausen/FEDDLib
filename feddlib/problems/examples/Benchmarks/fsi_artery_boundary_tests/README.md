# fsi_artery_boundary_tests

Boundary-condition benchmark using the 10 mm and 12 mm fluid/solid mesh pairs in
`meshes/artery_adaptations` and the original 8 mm pair in its `deprecated` folder.
The inlet uses a precomputed Laplace profile,
imported from HDF5 as in `curvedArteryFSI_Benchmark`.

## Build and run

From the workspace root, using an existing configured FEDDLib build:

```sh
cmake -S source/FEDDLib -B build
cmake --build build --target problems_fsi_artery_boundary_tests -j 4
cd build/feddlib/problems/examples/Benchmarks/fsi_artery_boundary_tests
mpirun -np 4 ./problems_fsi_artery_boundary_tests.exe --problemfile=parametersProblemFSI.xml
```

CMake copies the meshes, XML settings and five inlet profiles to the build
directory. P2 velocity profiles are provided for the original 8 mm fluid mesh,
fluid resolution 1 at 10 mm and 12 mm, and fluid resolutions 2 and 3 at 10 mm.
The default mesh pair is 12 mm;
to use the original 10 mm pair, change both `Mesh 1 Name` and `Mesh 2 Name`
to their `_10mm.mesh` variants. To use the resolution-2 fluid/solid pair, set:

```xml
<Parameter name="Mesh 1 Name" type="string" value="fsi_artery_fluid_res_2_Flags_4_5_6_9_10_15_10mm.mesh"/>
<Parameter name="Mesh 2 Name" type="string" value="fsi_solid_res_1_flags_6_7_8_9_10_11_12_13_14_15_10mm.mesh"/>
```

For the resolution-3 pair, set:

```xml
<Parameter name="Mesh 1 Name" type="string" value="fsi_artery_fluid_res_3_Flags_4_5_6_9_10_15_10mm.mesh"/>
<Parameter name="Mesh 2 Name" type="string" value="fsi_solid_res_3_flags_6_7_8_9_10_11_12_13_14_15_10mm.mesh"/>
```

For the original 8 mm pair, use the filenames without a length suffix:

```xml
<Parameter name="Mesh 1 Name" type="string" value="fsi_artery_fluid_res_1_Flags_4_5_6_9_10_15.mesh"/>
<Parameter name="Mesh 2 Name" type="string" value="fsi_solid_res_0_flags_6_7_8_9_10_11_12_13_14_15.mesh"/>
```

CMake copies these two meshes from `artery_adaptations/deprecated` into the build
directory. Their P2 inlet profile uses flag 4 at z = 0, directed along +z, like
the newer meshes; the outlet remains flag 5 at z = 8 mm.

The driver automatically selects
`laplace_parabolic_<fluid mesh basename>_<Discretization>.h5` from
`Mesh Partitioner/Mesh 1 Name` and `Parameter/Discretization`.
An optional `Parameter Fluid/Inflow Profile File` overrides this basename
(without the `.h5` extension). Its scalar field must be named `solution` and
use the same full fluid-node map and discretization as the benchmark.

This is a manually run benchmark; it has no CTest reference solution.

## Inlet profile and flow rate

`Parameter Fluid/BC Type` remains `parabolic_artery`. The spatial profile solves
`-Delta_Gamma p = 1` on inlet surface flag **4**, with `p = 0` on its rim.
It is stored with maximum value one and zero values away from the inlet.
All supplied fluid meshes use **z = 0** as the inlet, directed into the fluid
along **+z**, so the velocity callback sets `u = (0, 0, amplitude(t) * p)`.
The analytic ellipse parameters and `MaxVelocity` are no longer used.

The existing imported-vector boundary-condition mechanism integrates this
profile on the inlet mesh and scales its amplitude to
`Parameter Fluid/Flowrate` (cm³/s). A positive `Inflow Ramp Time` applies the
existing half-cosine ramp; zero disables the ramp. The same condition is
registered for the coupled FSI problem and its fluid component.
`Inflow Flag` must match `General/Flag Inlet Fluid` (4 for these profiles).
`Flowrate` also remains an input to the absorbing outlet model.

The inlet rim remains part of the FSI coupling, so its normal motion can add
to the prescribed profile's flux. The fluid wall is the moving solid interface.
The fluid outlet is **5** (`General/Flag Outlet Fluid`), interface **6**,
inlet/outlet rings **9/10**, solid ends **7/8**, and solid
support strips **13/14**.

Mesh loading uses `readAndPartition(15, "mm", true)`, converting millimetres
to centimetres. All physical parameters must use consistent units.

## Resistance outlet variants

Set these entries in `Parameter Fluid` to enable the relative model:

```xml
<Parameter name="Pressure Boundary Condition" type="string" value="Resistance"/>
<Parameter name="Resistance Type" type="string" value="Relative"/>
<Parameter name="Resistance" type="double" value="3.6660e-4"/>
<Parameter name="Distal Pressure" type="double" value="0.0"/>
<Parameter name="Unsteady Start" type="double" value="0.2"/>
```

`Resistance Type = Constant` is the default and uses the configured `Resistance`
as R0 throughout. `Relative` uses R0 before `Unsteady Start`; at and after that
time it uses
`alpha = 0.5 * (A_in(t)/A_in(t_u) + A_out(t)/A_out(t_u))` and `R(t) = R0/alpha^2`.
Here `t_u` is `Unsteady Start`. Both reference areas are captured once on the
first pressure evaluation at or after that time and retained across subsequent
timesteps and restarts. At capture, alpha is one and the resistance equals R0.
If the transition lies between time levels, capture occurs at the first time
level reaching it; no interpolation of areas is performed.
Both variants apply `p_out = Distal Pressure + R(t) * Q_used`; distal pressure
defaults to zero. `Reference fluid pressure` belongs to the absorbing models
and does not set the resistance model's distal pressure.

The existing explicit timestep update, flow extrapolation and `Average Flowrate`
option remain in use. For `Relative`, `Use Inflow = true` selects inlet flow
only before `Unsteady Start`; afterward outlet flow/history is used. For
`Constant`, `Use Inflow` retains its existing behavior throughout the run.
The pressure is evaluated at the start of a timestep from the retained flow
and current mesh; it has no additional Newton/Jacobian contribution.
The benchmark's default outlet remains `Absorbing Paper`.

## Regenerate the profiles

The generator uses FEDDLib's mesh maps, basis functions and HDF5 format. It
solves the inlet surface Poisson problem with conjugate gradients, checks the
residual and verifies an HDF5 import/export round trip. It does not run FSI.

```sh
cmake --build build --target problems_fsi_artery_inlet_profile -j 4
mesh="fsi_artery_fluid_res_1_Flags_4_5_6_9_10_15"
mpirun -np 1 build/feddlib/problems/examples/Benchmarks/fsi_artery_boundary_tests/problems_fsi_artery_inlet_profile.exe \
    --mesh="source/FEDDLib/meshes/artery_adaptations/deprecated/${mesh}.mesh" \
    --discretization=P2 \
    --inlet-flag=4 \
    --output="source/FEDDLib/feddlib/problems/examples/Benchmarks/fsi_artery_boundary_tests/laplace_parabolic_${mesh}_P2"
for length in 10 12; do
    mesh="fsi_artery_fluid_res_1_Flags_4_5_6_9_10_15_${length}mm"
    mpirun -np 1 build/feddlib/problems/examples/Benchmarks/fsi_artery_boundary_tests/problems_fsi_artery_inlet_profile.exe \
        --mesh="source/FEDDLib/meshes/artery_adaptations/${mesh}.mesh" \
        --discretization=P2 \
        --inlet-flag=4 \
        --output="source/FEDDLib/feddlib/problems/examples/Benchmarks/fsi_artery_boundary_tests/laplace_parabolic_${mesh}_P2"
done
for resolution in 2 3; do
    mesh="fsi_artery_fluid_res_${resolution}_Flags_4_5_6_9_10_15_10mm"
    mpirun -np 1 build/feddlib/problems/examples/Benchmarks/fsi_artery_boundary_tests/problems_fsi_artery_inlet_profile.exe \
        --mesh="source/FEDDLib/meshes/artery_adaptations/${mesh}.mesh" \
        --discretization=P2 \
        --inlet-flag=4 \
        --output="source/FEDDLib/feddlib/problems/examples/Benchmarks/fsi_artery_boundary_tests/laplace_parabolic_${mesh}_P2"
done
cmake --build build --target problems_fsi_artery_boundary_tests -j 4
```

When adding another geometry, generate a matching profile and add it to the
CMake copy list. Check the boundary labels, inlet direction in the velocity
callback, unit conversion and solid/ALE conditions for that geometry.
