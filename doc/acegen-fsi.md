# AceGen materials in FSI

Use a 3D P2 solid with the AceGen Interface enabled. In the FSI problem's
`Parameter Solid` list, select `Use AceGen Interface = true` and the desired
`Material model` (`SCI_NH`, `SCI_SMC_MLCK`, or
`SCI_SMC_Active_Growth_Reorientation`). `Structure Model` may be omitted; if
supplied, it must match `Material model`. The older configuration using
`Material model = Neo-Hooke`, `FSCI = true` (or `SCI = true`) and an explicit
`Structure Model` remains supported. FSI passes
the material and diffusion sublists to the nonlinear solid component.

Mechanics-only FSI initializes concentration to zero on the repeated solid map.
Call `NonLinElasticity::updateConcentration()` to prescribe a concentration;
the same field is used in assembly, post-processing, and history updates.
The solid commits its material history after an accepted FSI timestep.

## Passive restart and process markers

The element clock starts at the validated restart time and uses the continuation
increment selected by the time integrator. For example, restarting at `0.1`
with `dt = 0.001` evaluates the next material response at `0.101`.

Checkpoints for the MLCK and growth/reorientation materials contain a
`Material process state` sublist with three boolean markers: `Active response
took place`, `Growth took place`, and `Reorientation took place`. These stay
true once recorded, including after the corresponding process is switched off.
The markers are combined across all solid elements and MPI ranks and included
in recovery checkpoints as well. Native activation after `ActiveStartTime`
also counts as active history.
Entering a process interval marks it immediately, so checkpoints exactly at the
phase boundary are also excluded from passive restart.

Only passive checkpoints can currently be resumed. The passive reference fibre
state is initialized without resetting the physical clock. A saved true marker
causes a clear error before checkpoint fields are restored; full material-history
restoration remains separate work. Older checkpoints without these markers have
unknown material history and are rejected for these materials, including with
`Allow legacy restart = true`.

## Material output

For `SCI_SMC_Active_Growth_Reorientation`, enable the following in the FSI
problem's `Exporter` list:

```xml
<Parameter name="Export material fields" type="bool" value="true"/>
```

This writes `solid_material` ParaView output after convergence, using
`Export every X timesteps`. Fields include von Mises stress, circumferential,
axial and radial stress, energy, fibre directions, both stretches, and the
four `nC`/`nD` state variables. Projection sums the weighted nodal contributions
across MPI ranks and divides by the summed volume weights. Field columns use the original AceGen indices (10 for von Mises stress,
27/28 for the two stretches), with column zero supplying the volume weight. Unsupported materials return
empty data/name lists; explicitly requesting their export raises an error.

For a separate stress-only series with one value per solid element, set:

```xml
<Parameter name="Export element stresses" type="bool" value="true"/>
```

This writes `solid_stress_elements.xmf`/`.h5` with cell fields `vonMisesStress`,
`SCirc`, `SAxial`, and `SRadial`. It works independently of `Export material
fields` in FSI and SCI and follows the existing export cadence. Each value is
the current-volume-weighted Gauss-point average **before nodal recovery**;
von Mises stress is evaluated by AceGen at each Gauss point and then averaged.
The averaging uses a local copy of the material parameters with
`PostProcessingType = 1` and `PostProcessingPower = 0`, so it does not depend on
the configured nodal power or change material history. The existing P0 exporter
uses the four corner nodes of each tetrahedron for visualization.

Two additional, independent options are available in the same `Exporter` list:

```xml
<Parameter name="Export Gauss point stresses" type="bool" value="true"/>
<Parameter name="Export P1 stresses" type="bool" value="true"/>
```

Both default to false and support the growth/reorientation material in FSI and
SCI. `solid_stress_gauss.xmf`/`.h5` contains one point per integration point,
with the nine Cauchy stress components, `vonMisesStress`, zero-based `elementID`
and `gaussPointID`, and `currentVolumeWeight`. Coordinates are the reference
Gauss-point positions. The point-data vector `displacement` gives the interpolated
Gauss-point displacement; use ParaView's **Warp By Vector** with this field and
scale one to view the current positions. The points and stresses are evaluated
directly from AceGen task 6, without nodal averaging. The adapter computes
coordinates itself because the native task's coordinate-buffer indices overlap.
The external AceGen installation does not need modifications.

`solid_stress_p1.xmf`/`.h5` contains a continuous P1 L2 recovery of the same ten
stress quantities. It assembles the consistent mass matrix and right-hand side
from the current-volume Gauss weights and solves with diagonally preconditioned
CG. Von Mises stress is projected as its own scalar; it is not recomputed from
the recovered tensor. The output retains the original P2 tetrahedral geometry,
with every midside stress set to the mean of the two corner stresses. It therefore
avoids the different P2 corner/midside recovery weights and can be coloured as
point data in ParaView.

The recovery is continuous across shared nodes, including shared material
interfaces. It can smooth physical jumps and can overshoot Gauss-point extrema;
no clipping is applied. Keep the raw series for inspecting those values. Both
outputs follow the existing export cadence and leave material history unchanged.

## AceGen versions

The growth/reorientation adapter resolves material parameters by the names
reported by the installed library, since their ordering differs between
generated versions. Existing parameters keep their configured values/defaults.
`Gamma1Max` and `Gamma3Max` fall back to the older `Gamma1` and `Gamma3` values.
New parameters without an older equivalent must be supplied explicitly as
**double** entries in `Parameter Solid`; the error names any missing entry.
For the currently installed version these are `Gamma1Min`, `Gamma3Min`,
`ActiveBool`, `growthBool`, `reorientationBool`, `growthDirectionBool`,
`kFiberDistribution`, and `P`. Choose these values for your material model;
the regression test's values describe a passive test case.

The adapter enables nodal projection (`PostProcessingType = 1`) and accepts
`PostProcessingPower` in `Parameter Solid` (default 1).

Run `ctest -R aceGenStructure --output-on-failure` in the build directory for
the mechanics-only assembly, residual, post-processing and timestep regression.
