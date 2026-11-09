# fsi_artery_boundary_tests

Boundary-condition benchmark using the two fluid meshes in
`meshes/artery_adaptations`. The inlet uses a precomputed Laplace profile,
imported from HDF5 as in `curvedArteryFSI_Benchmark`.

## Build and run

From the workspace root, using an existing configured FEDDLib build:

```sh
cmake -S source/FEDDLib -B build
cmake --build build --target problems_fsi_artery_boundary_tests -j 4
cd build/feddlib/problems/examples/Benchmarks/fsi_artery_boundary_tests
mpirun -np 4 ./problems_fsi_artery_boundary_tests.exe --problemfile=parametersProblemFSI.xml
```

CMake copies the meshes, XML settings and four inlet profiles to the build
directory. Profiles are provided for fluid resolutions 1 and 2, each with P1
and P2 velocity discretization. The driver automatically selects
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
Both supplied meshes have an inlet at `z=0`, directed into the fluid along
**+z**, so the velocity callback sets `u = (0, 0, amplitude(t) * p)`.
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
The outlet, solid and ALE conditions are unchanged: fluid outlet **5**,
interface **6**, inlet/outlet rings **9/10**, solid ends **7/8**, and solid
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
`alpha = 0.5 * (A_in(t)/A_in(0) + A_out(t)/A_out(0))` and `R(t) = R0/alpha^2`.
The reference areas are taken from the initial mesh at time zero and retained
across restarts. They are not reset at the unsteady transition.
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
for resolution in 1 2; do
    mesh="fsi_artery_fluid_res_${resolution}_Flags_4_5_6_9_10_15"
    for fe in P1 P2; do
        mpirun -np 1 build/feddlib/problems/examples/Benchmarks/fsi_artery_boundary_tests/problems_fsi_artery_inlet_profile.exe \
            --mesh="source/FEDDLib/meshes/artery_adaptations/${mesh}.mesh" \
            --discretization="${fe}" \
            --output="source/FEDDLib/feddlib/problems/examples/Benchmarks/fsi_artery_boundary_tests/laplace_parabolic_${mesh}_${fe}"
    done
done
cmake --build build --target problems_fsi_artery_boundary_tests -j 4
```

When adding another geometry, generate a matching profile and add it to the
CMake copy list. Check the boundary labels, inlet direction in the velocity
callback, unit conversion and solid/ALE conditions for that geometry.
