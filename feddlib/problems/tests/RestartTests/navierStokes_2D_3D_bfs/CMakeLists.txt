TRIBITS_ADD_EXECUTABLE(
    navierStokes_2D_3D_bfs
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

# Each mesh/rank pair has its own clean directory. A four-rank run stops at 0.01;
# restart on four/six ranks from 0.005 to 0.01, then from the final checkpoint at
# 0.01 to 0.015 against an independent four-rank reference that stops at 0.015.
SET(reference_mpi_procs 4)
FOREACH(dimension 2 3)
    FOREACH(restart_mpi_procs 4 6)
        IF(NOT TPL_ENABLE_MPI OR reference_mpi_procs GREATER MPI_EXEC_MAX_NUMPROCS
            OR restart_mpi_procs GREATER MPI_EXEC_MAX_NUMPROCS)
            CONTINUE()
        ENDIF()
        SET(restart_test_name ${PACKAGE_NAME}_navierStokes_2D_3D_bfs_${dimension}D_${reference_mpi_procs}_to_${restart_mpi_procs})
        SET(restart_test_script ${CMAKE_CURRENT_BINARY_DIR}/${restart_test_name}.cmake)
        IF(dimension EQUAL 2)
            SET(restart_test_mesh BFS2d_1600.mesh)
        ELSE()
            SET(restart_test_mesh BFS3dCC.mesh)
        ENDIF()
        CONFIGURE_FILE(runRestartTest.cmake.in ${restart_test_script} @ONLY)
        ADD_TEST(NAME ${restart_test_name}
            COMMAND ${CMAKE_COMMAND}
                "-DEXECUTABLE=$<TARGET_FILE:${PACKAGE_NAME}_navierStokes_2D_3D_bfs>"
                -P ${restart_test_script}
        )
        SET_TESTS_PROPERTIES(${restart_test_name} PROPERTIES
            PROCESSORS ${restart_mpi_procs}
            LABELS ${PACKAGE_NAME}
            REQUIRED_FILES "$<TARGET_FILE:${PACKAGE_NAME}_navierStokes_2D_3D_bfs>"
        )
    ENDFOREACH()
ENDFOREACH()

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_navierStokes_2D_3D_bfs
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${navier_stokes_restart_parameter_files}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${navier_stokes_restart_parameter_files}
    EXEDEPS navierStokes_2D_3D_bfs
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_navierStokes_2D_3D_bfs
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES ${navier_stokes_restart_mesh_files}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${navier_stokes_restart_mesh_files}
    EXEDEPS navierStokes_2D_3D_bfs
)
