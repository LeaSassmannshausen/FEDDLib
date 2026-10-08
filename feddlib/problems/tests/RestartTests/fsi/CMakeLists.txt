TRIBITS_ADD_EXECUTABLE(fsi_artery_restart SOURCES main.cpp)

SET(fsi_restart_parameter_files
    parametersProblemFSI.xml parametersProblemFSI_restart.xml
    parametersPrecGE.xml parametersPrecFluidMono.xml parametersPrecStructure.xml
    parametersPrecGeometry.xml parametersSolverFSI.xml parametersSolverGeometry.xml
)
TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fsi_artery_restart
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${fsi_restart_parameter_files}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${fsi_restart_parameter_files}
    EXEDEPS fsi_artery_restart
)
TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_fsi_artery_restart
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES fsi_fluid_length_0_5_mm.mesh fsi_solid_length_0_5_mm.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES fsi_fluid_length_0_5_mm.mesh fsi_solid_length_0_5_mm.mesh
    EXEDEPS fsi_artery_restart
)

# Exercise all pressure models on the same 3D arterial segment. Each run has
# its own files and restores once before and once after the area transition.
IF(TPL_ENABLE_MPI AND MPI_EXEC_MAX_NUMPROCS GREATER_EQUAL 4)
    FOREACH(outlet_model Resistance Absorbing "Absorbing Paper")
        STRING(REPLACE " " "_" outlet_name "${outlet_model}")
        SET(outlet_test_name ${PACKAGE_NAME}_fsi_artery_restart_${outlet_name})
        SET(outlet_test_script ${CMAKE_CURRENT_BINARY_DIR}/${outlet_test_name}.cmake)
        CONFIGURE_FILE(runOutletRestartTest.cmake.in ${outlet_test_script} @ONLY)
        ADD_TEST(NAME ${outlet_test_name}
            COMMAND ${CMAKE_COMMAND}
                "-DEXECUTABLE=$<TARGET_FILE:${PACKAGE_NAME}_fsi_artery_restart>"
                -P ${outlet_test_script})
        SET_TESTS_PROPERTIES(${outlet_test_name} PROPERTIES
            PROCESSORS 4 LABELS ${PACKAGE_NAME} TIMEOUT 600
            REQUIRED_FILES "$<TARGET_FILE:${PACKAGE_NAME}_fsi_artery_restart>")
    ENDFOREACH()
ENDIF()
