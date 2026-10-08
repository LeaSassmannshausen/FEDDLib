TRIBITS_ADD_EXECUTABLE(
    unsteadyNavierStokes_restart
    SOURCES main.cpp
)

SET(navier_stokes_restart_parameter_files
    parametersProblem.xml
    parametersProblem2D.xml
    parametersProblem3D.xml
    parametersProblem_restart.xml
    parametersProblem_reference.xml
    parametersProblem_restart_final.xml
    parametersPrec.xml
    parametersSolver.xml
)
SET(navier_stokes_restart_mesh_files BFS2d_1600.mesh BFS3dCC.mesh)

# Each mesh/rank pair has its own clean directory. A four-rank run stops at 0.02;
# restart on four/six ranks from 0.01 to 0.02, then from the final checkpoint at
# 0.02 to 0.03 against an independent four-rank reference that stops at 0.03.
SET(reference_mpi_procs 4)
FOREACH(dimension 2 3)
    FOREACH(restart_mpi_procs 4 6)
        IF(NOT TPL_ENABLE_MPI OR reference_mpi_procs GREATER MPI_EXEC_MAX_NUMPROCS
            OR restart_mpi_procs GREATER MPI_EXEC_MAX_NUMPROCS)
            CONTINUE()
        ENDIF()
        SET(restart_test_name ${PACKAGE_NAME}_unsteadyNavierStokes_restart_${dimension}D_${reference_mpi_procs}_to_${restart_mpi_procs})
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
            PROCESSORS ${restart_mpi_procs}
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
