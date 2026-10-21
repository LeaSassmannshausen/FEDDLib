TRIBITS_ADD_EXECUTABLE(fsi_3D_arterialSegment SOURCES main.cpp)

SET(fsi_3D_arterialSegment_parameter_files
    parametersProblemFSI.xml parametersProblemFSI_restart.xml
    parametersPrecGE.xml parametersPrecFluidMono.xml parametersPrecStructure.xml
    parametersPrecGeometry.xml parametersSolverFSI.xml parametersSolverGeometry.xml
)
TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fsi_3D_arterialSegment
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${fsi_3D_arterialSegment_parameter_files}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${fsi_3D_arterialSegment_parameter_files}
    EXEDEPS fsi_3D_arterialSegment
)
TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_fsi_3D_arterialSegment
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES fsi_fluid_length_0_5_mm.mesh fsi_solid_length_0_5_mm.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES fsi_fluid_length_0_5_mm.mesh fsi_solid_length_0_5_mm.mesh
    EXEDEPS fsi_3D_arterialSegment
)

SET(arterial_reference_root
    ${CMAKE_CURRENT_SOURCE_DIR}/../../UnitTest/ReferenceSolutions/fsi_3D_arterialSegment)
FOREACH(reference_model Resistance Absorbing Absorbing_Paper)
    SET(reference_source ${arterial_reference_root}/${reference_model})
    FILE(GLOB reference_files RELATIVE ${reference_source}
        "${reference_source}/*.h5" "${reference_source}/*.xml")
    TRIBITS_COPY_FILES_TO_BINARY_DIR(references_fsi_3D_arterialSegment_${reference_model}
        SOURCE_DIR ${reference_source} SOURCE_FILES ${reference_files}
        DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions/${reference_model}
        DEST_FILES ${reference_files} EXEDEPS fsi_3D_arterialSegment)
ENDFOREACH()

# Exercise all pressure models from frozen checkpoints. Each run has its own
# files. Resistance has one save-all restart; absorbing models each have two
# restarts at two inlet-ramp times. Only Absorbing Paper captures an outlet area
# transition. Every simulation advances two steps.
IF(TPL_ENABLE_MPI AND MPI_EXEC_MAX_NUMPROCS GREATER_EQUAL 4)
    FOREACH(outlet_model Resistance Absorbing "Absorbing Paper")
        STRING(REPLACE " " "_" outlet_name "${outlet_model}")
        SET(outlet_test_name ${PACKAGE_NAME}_fsi_3D_arterialSegment_${outlet_name})
        SET(outlet_test_script ${CMAKE_CURRENT_BINARY_DIR}/${outlet_test_name}.cmake)
        CONFIGURE_FILE(runOutletRestartTest.cmake.in ${outlet_test_script} @ONLY)
        ADD_TEST(NAME ${outlet_test_name}
            COMMAND ${CMAKE_COMMAND}
                "-DEXECUTABLE=$<TARGET_FILE:${PACKAGE_NAME}_fsi_3D_arterialSegment>"
                -P ${outlet_test_script})
        SET_TESTS_PROPERTIES(${outlet_test_name} PROPERTIES
            PROCESSORS 4 LABELS ${PACKAGE_NAME} TIMEOUT 600
            REQUIRED_FILES "$<TARGET_FILE:${PACKAGE_NAME}_fsi_3D_arterialSegment>")
    ENDFOREACH()
ENDIF()

# Explicit maintenance operation; normal CTest cases never execute this script.
CONFIGURE_FILE(generateOutletReferences.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/generateOutletReferences.cmake @ONLY)
