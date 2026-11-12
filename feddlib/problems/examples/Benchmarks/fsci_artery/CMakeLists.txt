tribits_add_executable(fsci_artery SOURCES main.cpp)

set(fsci_parameters
    parametersProblemFSCI.xml parametersSolverFSCI.xml parametersSolverGeometry.xml
    parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecFluidBlock.xml
    parametersPrecStructure.xml parametersPrecGeometry.xml parametersPrecSCI.xml parametersPrecChem.xml
)
tribits_copy_files_to_binary_dir(data_fsciArtery
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${fsci_parameters}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${fsci_parameters}
    EXEDEPS fsci_artery
)

set(fsci_meshes
    fsi_artery_fluid_res_3_Flags_4_5_6_9_10_15_10mm.mesh
    fsi_solid_res_3_flags_6_7_8_9_10_11_12_13_14_15_10mm.mesh
)
tribits_copy_files_to_binary_dir(meshes_fsciArtery
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/artery_adaptations
    SOURCE_FILES ${fsci_meshes}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${fsci_meshes}
    EXEDEPS fsci_artery
)

set(fsci_legacyMeshes
    fsi_artery_fluid_res_1_Flags_4_5_6_9_10_15.mesh
    fsi_solid_res_0_flags_6_7_8_9_10_11_12_13_14_15.mesh
)
tribits_copy_files_to_binary_dir(legacyMeshes_fsciArtery
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/artery_adaptations/deprecated
    SOURCE_FILES ${fsci_legacyMeshes}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${fsci_legacyMeshes}
    EXEDEPS fsci_artery
)

set(fsci_profiles
    laplace_parabolic_fsi_artery_fluid_res_3_Flags_4_5_6_9_10_15_10mm_P2.h5
)
tribits_copy_files_to_binary_dir(profiles_fsciArtery
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../fsi_artery_boundary_tests
    SOURCE_FILES ${fsci_profiles}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${fsci_profiles}
    EXEDEPS fsci_artery
)
