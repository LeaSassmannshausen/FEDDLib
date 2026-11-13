# FSI/FSCI validation

This test runs FSI and implicit-chemistry FSCI independently on the simple arterial
segment used by `restartTests/fsi_3D_arterialSegment`. The coordinates remain in
centimetres (radius 0.09 cm, length 0.05 cm). Both cases use P2 velocity, solid and
ALE displacement, with P1 pressure, the same normalized circular inlet profile,
and a constant resistance outlet. Each run starts from zero and takes two steps
of 0.0025 s. The flow rate is 0.001 cm³/s with the restart fixture's cosine ramp.

Both problems use identical `SCI_SMC_Active_Growth_Reorientation` material
parameters. Activation, growth and reorientation are disabled. Diffusion is zero
in both the material and the transport problem; reaction and initial
concentration are zero, and there is no concentration inlet. The chemistry
equation remains part of the FSCI solve, and its concentration must remain zero.

At **each** accepted step, including the requested final comparison after two
steps, the test checks:

- Global L2 and Linf differences of `u_f`, `p`, `d_s`, `lambda` and `d_f`,
  importing by global ID. Each norm must satisfy
  `error <= absolute tolerance + relative tolerance * norm(FSI field)`.
  Output includes both field norms and their absolute and relative errors.
- The applied `pressureOutlet`, with bound `1e-10 + 1e-7 * |FSI pressureOutlet|`.
  This samples the existing explicit resistance load; it does not reassemble or
  advance the outlet state during comparison.
- Equal nonlinear counts and equal average/total linear iteration counts from
  the solver's iteration logs. `linearIterations` contains the average per
  Newton update; multiplying by `newtonIterations` gives the total.
- Two accepted steps at the requested times (`dt` and `2*dt`), matching clocks,
  zero chemistry, nonzero solver work, and nonzero final flow and outlet pressure.

The `Validation` sublist in `parameters.xml` controls the field/pressure
tolerances (absolute `1e-10`, relative `1e-7`) and the zero-concentration tolerance
(`1e-14`). Iteration counts must agree; these tolerances do not relax that check.

## Code structure

`main.cpp` follows the AssembleFE tests' setup/solve/compare pattern. Its final
three sections explicitly construct the two problems and compare their results:

1. **FSI PROBLEM** constructs `FSIValidationProblem fsi(...)`, attaches the
   mechanical boundary conditions, and advances two steps.
2. **FSCI PROBLEM** constructs `FSCIValidationProblem fsci(...)` with the same
   inputs plus the chemical domain, zero diffusion tensor and zero reaction.
   SCI's local solid/chemistry boundary setup is shown in this section.
3. **Comparison** checks both accepted steps, prints the errors and iteration
   counts, and returns failure if any check fails.

`FSIValidationProblem` and `FSCIValidationProblem` are aliases for small subclasses
of the production FSI/FSCI classes. Their only override records accepted
solutions before the timestep history is committed. They do not change assembly
or the solution algorithm. Each recorded vector is a deep copy.

The common helpers above `main()` handle component parameter lists, independent
P1/P2 domains, the temporal block definition, mechanical boundary conditions,
time integration and comparisons. Comments explain the mesh flags, GE/GI block
ordering, inlet flux normalization, pressure sampling and iteration log format.

All physical and solver inputs are shared. The formulation flags, the extra
chemistry block and FaCSI/FaCSCI preconditioning differ. Component solver inputs
and meshes are copied from the existing arterial restart fixture by CMake.
No production classes are modified by this test.

Build and run from the build directory:

```sh
cmake --build . --target problems_FSI_FSCI_Validation_Test -j4
ctest --output-on-failure -R FSI_FSCI_Validation_Test
```

CTest registers one- and four-rank cases when available, with AceGen Interface
enabled. The default uses explicit geometry; the executable also accepts
`--geometry-implicit`. Iteration logs use formulation/geometry/rank suffixes so
the registered MPI runs can execute concurrently.

As in the AssembleFE drivers, input files can be selected with `--problemfile`,
`--solverfile`, `--precfileFluid`, `--precfileStructure`, `--precfileGeometry` and
`--solverfileGeometry`; `--help` lists them. Custom physical inputs must retain
the passive, zero-diffusion setup and exactly two timesteps from zero.
