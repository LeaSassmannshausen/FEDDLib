tribits_add_executable(sci_artery SOURCES main.cpp)

set(sci_artery_parameters
    parametersProblemSCI.xml parametersSolverSCI.xml
    parametersPrecStructure.xml parametersPrecChem.xml parametersPrecCE.xml
)
tribits_copy_files_to_binary_dir(data_sciArtery
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${sci_artery_parameters}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${sci_artery_parameters}
    EXEDEPS sci_artery
)
set(sci_artery_meshes
    fsi_solid_res_0_flags_6_7_8_9_10_11_12_13_14_15_10mm.mesh
    fsi_solid_res_0_flags_6_7_8_9_10_11_12_13_14_15_12mm.mesh
)
tribits_copy_files_to_binary_dir(meshes_sciArtery
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/artery_adaptations
    SOURCE_FILES ${sci_artery_meshes}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${sci_artery_meshes}
    EXEDEPS sci_artery
)
