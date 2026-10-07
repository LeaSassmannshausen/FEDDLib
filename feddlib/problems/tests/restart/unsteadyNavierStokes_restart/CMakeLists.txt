TRIBITS_ADD_EXECUTABLE(
    unsteadyNavierStokes_restart
    SOURCES main.cpp
)

SET(navier_stokes_restart_parameter_files
    parametersProblem.xml
    parametersProblem2D.xml
    parametersProblem3D.xml
    parametersProblem_restart.xml
    parametersPrec.xml
    parametersSolver.xml
)
SET(navier_stokes_restart_mesh_files BFS2d_1600.mesh BFS3dCC.mesh)

# Each mesh/rank combination has its own clean directory. Phase 1 generates the
# checkpoints at 0.01 and 0.02; phase 2 resumes at 0.01 and compares both fields at 0.02.
FOREACH(dimension 2 3)
    FOREACH(mpi_procs 4 6)
        IF(NOT TPL_ENABLE_MPI OR mpi_procs GREATER MPI_EXEC_MAX_NUMPROCS)
            CONTINUE()
        ENDIF()
        SET(restart_test_name ${PACKAGE_NAME}_unsteadyNavierStokes_restart_${dimension}D_${mpi_procs})
        SET(restart_test_script ${CMAKE_CURRENT_BINARY_DIR}/${restart_test_name}.cmake)
        IF(dimension EQUAL 2)
            SET(restart_test_mesh BFS2d_1600.mesh)
        ELSE()
            SET(restart_test_mesh BFS3dCC.mesh)
        ENDIF()
        CONFIGURE_FILE(runRestartTest.cmake.in ${restart_test_script} @ONLY)
        ADD_TEST(NAME ${restart_test_name}
            COMMAND ${CMAKE_COMMAND}
                "-DEXECUTABLE=$<TARGET_FILE:${PACKAGE_NAME}_unsteadyNavierStokes_restart>"
                -P ${restart_test_script}
        )
        SET_TESTS_PROPERTIES(${restart_test_name} PROPERTIES
            PROCESSORS ${mpi_procs}
            LABELS ${PACKAGE_NAME}
            REQUIRED_FILES "$<TARGET_FILE:${PACKAGE_NAME}_unsteadyNavierStokes_restart>"
        )
    ENDFOREACH()
ENDFOREACH()

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadyNavierStokes_restart
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${navier_stokes_restart_parameter_files}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${navier_stokes_restart_parameter_files}
    EXEDEPS unsteadyNavierStokes_restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadyNavierStokes_restart
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES ${navier_stokes_restart_mesh_files}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${navier_stokes_restart_mesh_files}
    EXEDEPS unsteadyNavierStokes_restart
)
