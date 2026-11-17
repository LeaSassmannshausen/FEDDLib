IF(TPL_ENABLE_AceGENInterface)
    TRIBITS_ADD_EXECUTABLE(FSI_FSCI_Validation_Test SOURCES main.cpp)
    TRIBITS_COPY_FILES_TO_BINARY_DIR(data_FSI_FSCI_Validation_Test
        SOURCE_FILES parameters.xml parametersSolver.xml
        EXEDEPS FSI_FSCI_Validation_Test)
    # Share the arterial restart fixture, including its component solvers.
    TRIBITS_COPY_FILES_TO_BINARY_DIR(solvers_FSI_FSCI_Validation_Test
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../restartTests/fsi_3D_arterialSegment
        SOURCE_FILES parametersPrecFluidMono.xml parametersPrecStructure.xml
            parametersPrecGeometry.xml parametersSolverGeometry.xml
        EXEDEPS FSI_FSCI_Validation_Test)
    TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_FSI_FSCI_Validation_Test
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
        SOURCE_FILES fsi_fluid_length_0_5_mm.mesh fsi_solid_length_0_5_mm.mesh
        EXEDEPS FSI_FSCI_Validation_Test)
    FOREACH(ranks 1 4)
        # GI adds ALE work at each Newton update; allow for single-rank LU solves.
        TRIBITS_ADD_TEST(FSI_FSCI_Validation_Test
            NAME tests_FSI_FSCI_Validation_Test_GI
            ARGS "--geometry-implicit"
            NUM_MPI_PROCS ${ranks}
            TIMEOUT 900)
    ENDFOREACH()
    FOREACH(geometry implicit explicit)
        TRIBITS_ADD_TEST(FSI_FSCI_Validation_Test
            NAME tests_FSCI_passive_restart_${geometry}
            ARGS "--passive-restart --geometry-${geometry}"
            NUM_MPI_PROCS 4
            TIMEOUT 900)
    ENDFOREACH()
    feddlib_set_test_runtimes(DEFAULT High)
ENDIF()
