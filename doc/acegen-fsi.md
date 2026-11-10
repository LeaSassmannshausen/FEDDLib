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
