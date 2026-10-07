# Checkpointing and restart

Time-dependent problems can write HDF5 checkpoints and resume from a saved
solution. This branch provides the shared solution/BDF/Newmark checkpoint code
and FSI restart support. SCI/FSCI solvers and material-history extensions are
separate features.

Configure an uninterrupted run with these entries in `Timestepping Parameter`:

```xml
<Parameter name="Checkpointing" type="bool" value="true"/>
<Parameter name="Checkpoint directory" type="string" value="restartCheckpoints"/>
<Parameter name="Number Checkpoints" type="int" value="2"/>
<ParameterList name="Checkpoints">
    <Parameter name="1" type="double" value="0.01"/>
    <Parameter name="2" type="double" value="0.02"/>
</ParameterList>
```

Resume using the same mesh, discretization, variable names, time integration
scheme, and time step size, with a later `Final time`:

```xml
<Parameter name="Checkpointing" type="bool" value="false"/>
<Parameter name="Restart" type="bool" value="true"/>
<Parameter name="Restart directory" type="string" value="restartCheckpoints"/>
<Parameter name="Time step" type="double" value="0.01"/>
```

`Time step` is the physical restart time; `dt` is the time increment. Both
directory parameters default to the working directory. Checkpointing and restart
are disabled by default. Use separate output and input directories if a resumed
run also writes checkpoints, since creating an output file truncates it.

Checkpoint times should coincide with time steps. For BDF/FSI, checkpoints are
written when the next step begins; let the uninterrupted run advance one step
past a checkpoint you need. A multistep restart needs the saved previous steps
as well as the solution at the restart time. The tests restart after several
steps, so that this history is available.

The checkpoint files include `Solution<variable>.h5`, Newmark displacement,
velocity and acceleration, and, for FSI, moving-mesh mass products and geometry.
Dataset names are physical times formatted by `std::to_string`. Transfer the
complete checkpoint directory, rather than just the displacement or velocity
file. The FSI test covers explicit geometry with the Turek benchmark; it does
not validate restart of stateful outlet boundary conditions.

## Validation

With tests enabled and the required Trilinos solvers available:

```sh
cmake --build build --target problems_unsteadyNavierStokes_Restart problems_fsi_restart
ctest --test-dir build --output-on-failure -R 'problems_(unsteadyNavierStokes_Restart|fsi_restart)'
```

- `unsteadyNavierStokes_Restart` reads the checked-in HDF5 fixture and compares
  the continued velocity solution at the final time, on four and six MPI ranks.
- `fsi_restart` runs an uninterrupted simulation and a restarted simulation on
  two MPI ranks, comparing fluid velocity, pressure, and solid displacement at
  the same final time with relative tolerance `1e-8`. Both phases load the case
  from `parametersProblemFSI.xml`; `parametersProblemFSI_restart.xml` contains
  only the overrides to resume at `0.01` and stop at `0.02`. The test disables
  visualization and benchmark exports and reports the relative error for each
  of the three fields.

The restart tests use upstream meshes and do not require the new SCI meshes or
the AceGen Interface2 material models.
