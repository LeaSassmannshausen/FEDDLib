# Test runtime tiers

The tiers below use the user-provided laptop CTest results from 2026-10-21.
A timing covers the complete CTest case, including all MPI launches, setup,
reference production and validation phases. The categories remain fixed when
running on another machine; a later run does not automatically reclassify tests.
They describe expected wall-clock time, not enforced time limits.

| Test tier | Baseline runtime | Cases in this baseline |
|---|---|---|
| Low | Less than 10 s | 92 |
| Medium | At least 10 s and less than 30 s | 26 |
| High | At least 30 s | 6 |

## Configure-time selection

Enable tests as usual and pass `-D FEDDLib_TEST_RUNTIME:STRING=Low`, `Medium`,
or `High` in the configure script. Selection is cumulative:

| Configure value | Enabled tiers | Baseline cases enabled |
|---|---|---|
| `Low` | Low | 92 |
| `Medium` | Low and Medium | 118 |
| `High` (default) | Low, Medium and High | 124 |

These counts assume all baseline packages and up to nine MPI ranks are enabled.
Other package/MPI settings can register fewer tests. Invalid tier values are
rejected during configuration.

The supplied log names the general FSI example `problems_fsi_MPI_2` (20.80 s).
Its current source registration uses four ranks (`problems_fsi_MPI_4`); it keeps
the Medium estimate from that case. The original two-rank timing is retained
in the baseline table below and should be updated after measuring the new run.

For an existing workspace with `source/FEDDLib` and `build` directories:

```sh
cmake -S source/FEDDLib -B build -D FEDDLib_TEST_RUNTIME:STRING=Medium
ctest --test-dir build --output-on-failure
```

Reconfigure with `High` to enable all runtime tiers again. CTest still lists
excluded tests, showing them as disabled; their executables are not launched.
The option controls execution, while the test build targets remain available.
Existing package, MPI-rank, and other test-disable settings still apply.

The sample configure script exposes the same setting as `TEST_RUNTIME`, which
can be set in the script or its environment; the default is `High`.

## Labels and test maintenance

Every classified case additionally has exactly one of `RuntimeLow`,
`RuntimeMedium`, or `RuntimeHigh` appended to its existing CTest labels.
For example, with all tiers enabled:

```sh
ctest --test-dir build -L '^RuntimeLow$' --output-on-failure
ctest --test-dir build -L '^Runtime(Low|Medium)$' --output-on-failure
```

Labels provide further filtering within the configured selection; they do not
re-enable tests excluded by the configure-time setting.

Each test-registration `CMakeLists.txt` calls `feddlib_set_test_runtimes()` after
registering its tests. `DEFAULT` sets the tier for that directory, while `LOW`,
`MEDIUM`, and `HIGH` lists override individual full CTest names. Rank variants
can consequently have different tiers. For example:

```cmake
feddlib_set_test_runtimes(
    DEFAULT Low
    MEDIUM ${PACKAGE_NAME}_unsteadyNonLinElasticity_restart_2D_4_to_6
)
```

New test directories must also call this helper; configuration checks that every
registered test was classified. Without an explicit `DEFAULT`, the helper uses
`High`. Override names for tests not registered in the current package/MPI
configuration are ignored. Update the baseline below when reviewing new timings.

## Laptop baseline

| CTest case | Runtime (s) | Tier |
|---|---:|---|
| `core_domain_MPI_4` | 1.47 | Low |
| `core_elements_MPI_4` | 1.36 | Low |
| `core_fe_MPI_4` | 0.72 | Low |
| `core_quadRules_MPI_1` | 1.33 | Low |
| `core_testQuadRules_MPI_1` | 1.08 | Low |
| `core_surfaceIntegral_MPI_2` | 1.40 | Low |
| `core_neumannBC_MPI_4` | 1.44 | Low |
| `core_writeRead_MPI_4` | 1.42 | Low |
| `core_parameterList_MPI_1` | 1.07 | Low |
| `core_AABBTree_MPI_1` | 1.34 | Low |
| `core_mesh_AABBTree_MPI_1` | 1.39 | Low |
| `core_mesh_structured_MPI_4` | 1.39 | Low |
| `core_mesh_structured_3D_5_element_MPI_8` | 1.42 | Low |
| `core_mesh_elementFlags_MPI_2` | 2.77 | Low |
| `core_mesh_nodeFlagsP2_MPI_4` | 1.43 | Low |
| `core_mesh_surfaceOrientation_MPI_4` | 1.63 | Low |
| `core_mesh_unstructured_preProcessing_3D_cube_h_1_MPI_4` | 1.41 | Low |
| `core_mesh_unstructured_preProcessing_3D_fluid_mesh_MPI_4` | 2.23 | Low |
| `core_mesh_unstructured_preProcessing_2D_rectangle_MPI_4` | 1.49 | Low |
| `core_mesh_unstructured_writeMesh_3D_P2_MPI_3` | 1.37 | Low |
| `core_mesh_unstructured_MPI_2` | 1.45 | Low |
| `core_mesh_unstructured_writeMesh_2D_P2_MPI_2` | 1.42 | Low |
| `core_mesh_unstructured_writeMesh_2D_P1_MPI_4` | 1.44 | Low |
| `core_mesh_unstructured_writeMesh_3D_P1_MPI_5` | 1.46 | Low |
| `core_meshes_interface_MPI_2` | 1.67 | Low |
| `core_meshes_custom_partition_MPI_2` | 1.73 | Low |
| `core_mesh_flag_write_MPI_2` | 2.98 | Low |
| `core_map_tpetra_MPI_4` | 0.72 | Low |
| `core_matrix_MPI_4` | 1.36 | Low |
| `core_blockMatrixPtrAccess_MPI_2` | 1.36 | Low |
| `core_matrix_assembly_MPI_2` | 1.36 | Low |
| `core_multivector_MPI_4` | 1.34 | Low |
| `core_fe_test_laplace_MPI_4` | 1.35 | Low |
| `core_fe_test_navierstokes_MPI_4` | 1.42 | Low |
| `core_fe_test_linElas_MPI_1` | 1.37 | Low |
| `core_fe_test_nonlinElas_MPI_4` | 1.40 | Low |
| `core_interface_test_MPI_1` | 1.05 | Low |
| `problems_facsi_MPI_2` | 1.33 | Low |
| `problems_steadyNavierStokesAssFE_2D_unstructured_P2_P1_MPI_4` | 4.90 | Low |
| `problems_steadyNavierStokesAssFE_2D_structured_P1_P1_MPI_9` | 3.26 | Low |
| `problems_steadyNavierStokesAssFE_2D_structured_P2_P1_MPI_9` | 5.06 | Low |
| `problems_steadyNavierStokesAssFE_3D_structured_P1_P1_MPI_9` | 4.70 | Low |
| `problems_steadyNavierStokesAssFE_3D_structured_P2_P1_MPI_9` | 7.99 | Low |
| `problems_steadyNavierStokesAssFE_2D_structured_P2_P1_Newton_MPI_9` | 4.53 | Low |
| `problems_steadyNavierStokesAssFE_3D_structured_P2_P1_Newton_MPI_9` | 8.15 | Low |
| `problems_unsteadyNavierStokesAssFE_MPI_4` | 6.96 | Low |
| `problems_steadyGeneralizedNewtonianAssFE_MPI_4` | 20.62 | Medium |
| `problems_unsteadyLinElasAssFE_MPI_2` | 1.69 | Low |
| `problems_steadyLinElasAssFE_MPI_4` | 1.60 | Low |
| `problems_steadyNonLinElasAssFE_MPI_4` | 2.44 | Low |
| `problems_timeStepping_MPI_2` | 1.33 | Low |
| `problems_checkpointMetadata_MPI_2` | 3.80 | Low |
| `problems_recoveryCheckpoint_MPI_2` | 1.48 | Low |
| `problems_outputHistory_MPI_2` | 3.58 | Low |
| `problems_fsi_2D_P2P1_MPI_4` | 14.04 | Medium |
| `problems_unsteadyNonLinElasticity_2D_P1_4cores_MPI_4` | 8.99 | Low |
| `problems_unsteadyNonLinElasticity_2D_P1_6cores_MPI_6` | 10.39 | Medium |
| `problems_laplace_2D_P1_MPI_4` | 1.09 | Low |
| `problems_laplace_2D_P2_MPI_4` | 1.39 | Low |
| `problems_laplace_3D_P1_MPI_8` | 1.67 | Low |
| `problems_laplace_3D_P2_MPI_8` | 21.68 | Medium |
| `problems_linElas_2D_P1_MPI_4` | 0.76 | Low |
| `problems_linElas_2D_P2_MPI_4` | 1.38 | Low |
| `problems_linElas_3D_P1_MPI_8` | 2.11 | Low |
| `problems_linElas_3D_P2_MPI_8` | 43.33 | High |
| `problems_stokes_2D_P1_MPI_4` | 0.93 | Low |
| `problems_stokes_2D_P2_MPI_4` | 1.39 | Low |
| `problems_stokes_3D_P1_MPI_8` | 1.82 | Low |
| `problems_stokes_3D_P2_MPI_8` | 6.98 | Low |
| `problems_GNF_2D_P2_MPI_4` | 11.37 | Medium |
| `problems_nonLinElasticity_2D_P1_MPI_4` | 0.95 | Low |
| `problems_nonLinElasticity_3D_P1_MPI_8` | 19.58 | Medium |
| `problems_nonLinLaplace_2D_P1_MPI_4` | 0.95 | Low |
| `problems_nonLinLaplace_2D_P1_MPI_8` | 18.67 | Medium |
| `problems_navierStokes_2D_P1_Newton_MPI_4` | 0.85 | Low |
| `problems_navierStokes_2D_P2_Newton_MPI_4` | 1.57 | Low |
| `problems_navierStokes_3D_P1_Newton_MPI_8` | 2.94 | Low |
| `problems_navierStokes_3D_P2_Newton_MPI_8` | 45.74 | High |
| `problems_navierStokes_2D_P1_Monolithic_MPI_9` | 1.13 | Low |
| `problems_navierStokes_2D_P2_Block_MPI_9` | 2.68 | Low |
| `problems_navierStokes_2D_P2_PCD_MPI_9` | 2.97 | Low |
| `problems_navierStokes_2D_P2_Teko_MPI_9` | 2.07 | Low |
| `problems_fsi_3D_arterialSegment_Resistance` | 21.98 | Medium |
| `problems_fsi_3D_arterialSegment_Absorbing` | 21.04 | Medium |
| `problems_fsi_3D_arterialSegment_Absorbing_Paper` | 20.83 | Medium |
| `problems_fsi_2D_turek` | 23.00 | Medium |
| `problems_navierStokes_2D_3D_bfs_2D_4_to_4` | 27.88 | Medium |
| `problems_navierStokes_2D_3D_bfs_2D_4_to_6` | 18.58 | Medium |
| `problems_navierStokes_2D_3D_bfs_3D_4_to_4` | 28.80 | Medium |
| `problems_navierStokes_2D_3D_bfs_3D_4_to_6` | 27.44 | Medium |
| `problems_unsteadyNonLinElasticity_restart_2D_4_to_4` | 9.87 | Low |
| `problems_unsteadyNonLinElasticity_restart_2D_4_to_6` | 10.04 | Medium |
| `problems_fsi_2D_turek_recovery` | 62.10 | High |
| `problems_navierStokes_2D_3D_bfs_recovery` | 25.68 | Medium |
| `problems_navierStokes_2D_3D_initialSolution_2D_4_to_4` | 15.38 | Medium |
| `problems_navierStokes_2D_3D_initialSolution_2D_4_to_6` | 12.64 | Medium |
| `problems_navierStokes_2D_3D_initialSolution_3D_4_to_4` | 17.18 | Medium |
| `problems_navierStokes_2D_3D_initialSolution_3D_4_to_6` | 18.41 | Medium |
| `problems_navierStokes_variableTimeStep_restart` | 16.83 | Medium |
| `problems_fsi_variableTimeStep_restart` | 49.88 | High |
| `problems_newmark_variableTimeStep_restart` | 20.14 | Medium |
| `problems_artery_variableTimeStep_restart` | 55.21 | High |
| `problems_navierStokes_bfs_performance_2D_MPI_9` | 5.33 | Low |
| `problems_navierStokes_bfs_performance_3D_MPI_9` | 10.77 | Medium |
| `problems_navierStokes_bfs_performance_parameterFile_MPI_9` | 3.49 | Low |
| `problems_fsi_MPI_2` | 20.80 | Medium |
| `problems_geometry_MPI_4` | 2.18 | Low |
| `problems_laplace_MPI_4` | 1.70 | Low |
| `problems_steadyLinElas_MPI_4` | 1.72 | Low |
| `problems_steadyLinElas_Perf_MPI_8` | 2.12 | Low |
| `problems_unsteadyReactionDiffusion_MPI_4` | 2.27 | Low |
| `problems_steadyNavierStokes_MPI_4` | 2.01 | Low |
| `problems_steadyGeneralizedNewtonianAssFE_PowerLaw_MPI_4` | 43.81 | High |
| `problems_stokes_MPI_4` | 1.28 | Low |
| `problems_unsteadyLinElas_MPI_2` | 1.79 | Low |
| `problems_unsteadyNavierStokes_MPI_2` | 5.85 | Low |
| `problems_unsteadyNonLinElasticity_MPI_4` | 2.17 | Low |
| `amr_laplaceAdaptive2D_MPI_4` | 4.13 | Low |
| `amr_laplaceAdaptive3D_MPI_4` | 6.85 | Low |
| `amr_stokesAdaptive_MPI_2` | 8.17 | Low |
| `amr_steadyNavierStokesAdaptive_MPI_4` | 3.66 | Low |
| `amr_refinementFlagsTestUniform_MPI_4` | 24.73 | Medium |
| `amr_refinementFlagsTestArea_MPI_4` | 1.64 | Low |
| `amr_unsteadyReactionDiffusion_MPI_4` | 2.28 | Low |
