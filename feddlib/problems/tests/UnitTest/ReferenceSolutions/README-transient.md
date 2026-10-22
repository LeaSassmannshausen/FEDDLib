# FSI and unsteady nonlinear elasticity references

These HDF5 files are frozen numerical regression baselines, generated from an
uninterrupted simulation on four MPI ranks. They use the existing unit-test
format: one `solution` dataset per file. Normal tests only read these files;
missing or changed references cannot be replaced automatically by the test.

| Case | Comparison time | Fields | CTest ranks |
| --- | --- | --- | --- |
| `fsiUnitTest` | 0.01 | Fluid velocity (P2), pressure (P1), solid displacement (P2) | 4 |
| `unsteadyNonLinElasticityUnitTest` | 0.02 | Solid displacement, velocity and acceleration (P1) | 4 and 6 |

The FSI case uses the Turek fluid/solid `h004` meshes, explicit geometry and
FaCSI. Its independent inputs are the six `UnitTest/parameters_FSI_*.xml` files.
It takes four steps with `dt = 0.0025`, stops at `0.01` and compares the primary
fields. The nonlinear residual tolerance and coupled linear solver tolerance
are both `1e-8`; the geometry linear solver tolerance is `1e-10`.

The elasticity case uses `square_h02.mesh`, P1 elements, Saint Venant-Kirchhoff material
with Poisson ratio `0.48`,
a clamped left edge and a constant right-edge surface load, as in
`restartTests/unsteadyNonLinElasticity_restart`. Newmark parameters are `beta = 0.25`,
`gamma = 0.5`, `dt = 0.0025`. It runs to `0.0225`: the existing Newmark loop
finalizes derivatives at the start of the next step. The compared displacement
history, velocity and acceleration buffers therefore all represent `0.02`.
Changing this update convention requires adjusting which buffers are compared.

CMake copies the independent FSI inputs and the elasticity restart inputs into
the unit-test build directory. The
elasticity input names receive an `UnsteadyNonLinElasticity` suffix to avoid
overwriting other unit tests' solver files. Checkpointing, restart and save-all
output are disabled in the FSI input file and by the elasticity driver.
No restart simulation is performed.

Every field must have a finite, nonzero norm. Each comparison requires both
relative l2 error at most `1e-12` and absolute infinity error at most
`1e-11 + 1e-12 * ||reference||_infinity`.
The four-rank elasticity reference is also used for the six-rank test.

Initial validation: all three registered tests pass. Four-rank comparisons
report zero error; the largest six-rank relative error is approximately
`2.52e-14` (acceleration), with infinity error `2.81e-12`. Scaling each field
in temporary reference copies by `1.1` makes both drivers reject all three
comparisons and exit with failure. The source reference files remain unchanged.

From the repository workspace containing `source/FEDDLib` and `build`:

```sh
cmake --build build --target problems_fsiUnitTest problems_unsteadyNonLinElasticityUnitTest -j 2
ctest --test-dir build -V -R 'problems_(fsi_2D_P2P1|unsteadyNonLinElasticity_2D_P1_)'
```

To deliberately regenerate baselines after reviewing a numerical change, run
the executables with `--write-reference` and an explicit output directory.
For example, from `build/feddlib/problems/tests/UnitTest`:

```sh
mkdir -p /tmp/transient-unit-new-references
mpirun -np 4 ./problems_fsiUnitTest.exe --write-reference --reference-directory=/tmp/transient-unit-new-references
mpirun -np 4 ./problems_unsteadyNonLinElasticityUnitTest.exe --write-reference --reference-directory=/tmp/transient-unit-new-references
```

Review the generated files before replacing the six matching files in this
source directory, then rebuild to copy them into the build tree. The `4cores`
suffix identifies the generating run, including when a file is read on six
ranks. These baselines detect changes in the computed solutions; they are not
an independent analytic verification of the physical model.
