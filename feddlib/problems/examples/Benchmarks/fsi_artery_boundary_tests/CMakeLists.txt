tribits_add_executable(
    fsi_artery_boundary_tests
    SOURCES main.cpp
)

set(boundary_test_parameters
    parametersPrecGE.xml
    parametersPrecGI.xml
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

tribits_copy_files_to_binary_dir(meshes_fsiArteryBoundaryTests
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES arterie_with_plaque_initma_media_adventitia_fluid_fine.mesh arterie_with_plaque_initma_media_adventitia_solid_fine.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES arterie_with_plaque_initma_media_adventitia_fluid_fine.mesh arterie_with_plaque_initma_media_adventitia_solid_fine.mesh
    EXEDEPS fsi_artery_boundary_tests
)

set(artery_adaptation_meshes
    fsi_artery_fluid_res_1_Flags_4_5_6_9_10_15.mesh
    fsi_artery_fluid_res_2_Flags_4_5_6_9_10_15.mesh
    fsi_solid_res_1_flags_6_7_8_9_10_11_12_13_14_15.mesh
    fsi_solid_res_2_flags_6_7_8_9_10_11_12_13_14_15.mesh
    fsi_solid_res_3_flags_6_7_8_9_10_11_12_13_14_15.mesh
)

tribits_copy_files_to_binary_dir(adaptedMeshes_fsiArteryBoundaryTests
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/artery_adaptations
    SOURCE_FILES ${artery_adaptation_meshes}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${artery_adaptation_meshes}
    EXEDEPS fsi_artery_boundary_tests
)
