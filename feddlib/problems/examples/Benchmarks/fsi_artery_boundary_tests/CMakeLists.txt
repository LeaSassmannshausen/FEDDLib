tribits_add_executable(
    fsi_artery_boundary_tests
    SOURCES main.cpp
)

tribits_add_executable(
    fsi_artery_inlet_profile
    SOURCES inlet_profile.cpp
)

set(boundary_test_parameters
    parametersPrecGeometry.xml
    parametersProblemFSI.xml
    parametersSolverFSI.xml
    parametersSolverGeometry.xml
    parametersPrecFluidMono.xml
    parametersPrecFluidTeko.xml
    parametersPrecFluidBlock.xml
    parametersPrecStructure.xml
)

tribits_copy_files_to_binary_dir(data_fsiArteryBoundaryTests
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${boundary_test_parameters}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${boundary_test_parameters}
    EXEDEPS fsi_artery_boundary_tests
)

set(artery_adaptation_meshes
    fsi_artery_fluid_res_1_Flags_4_5_6_9_10_15_10mm.mesh
    fsi_artery_fluid_res_1_Flags_4_5_6_9_10_15_12mm.mesh
    fsi_artery_fluid_res_2_Flags_4_5_6_9_10_15_10mm.mesh
    fsi_solid_res_0_flags_6_7_8_9_10_11_12_13_14_15_10mm.mesh
    fsi_solid_res_0_flags_6_7_8_9_10_11_12_13_14_15_12mm.mesh
    fsi_solid_res_1_flags_6_7_8_9_10_11_12_13_14_15_10mm.mesh
)

tribits_copy_files_to_binary_dir(adaptedMeshes_fsiArteryBoundaryTests
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/artery_adaptations
    SOURCE_FILES ${artery_adaptation_meshes}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${artery_adaptation_meshes}
    EXEDEPS fsi_artery_boundary_tests
)

set(artery_inlet_profiles
    laplace_parabolic_fsi_artery_fluid_res_1_Flags_4_5_6_9_10_15_10mm_P2.h5
    laplace_parabolic_fsi_artery_fluid_res_1_Flags_4_5_6_9_10_15_12mm_P2.h5
    laplace_parabolic_fsi_artery_fluid_res_2_Flags_4_5_6_9_10_15_10mm_P2.h5
)

tribits_copy_files_to_binary_dir(inletProfiles_fsiArteryBoundaryTests
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${artery_inlet_profiles}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${artery_inlet_profiles}
    EXEDEPS fsi_artery_boundary_tests
)
