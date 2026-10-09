TRIBITS_ADD_EXECUTABLE(navierStokes_2D_3D_initialSolution SOURCES main.cpp)

# Share the BFS meshes and solver settings with the existing restart cases.
SET(bfs_case_directory ${CMAKE_CURRENT_SOURCE_DIR}/../navierStokes_2D_3D_bfs)
TRIBITS_COPY_FILES_TO_BINARY_DIR(data_navierStokes_2D_3D_initialSolution
    SOURCE_DIR ${bfs_case_directory}
    SOURCE_FILES parametersProblem2D.xml parametersProblem3D.xml parametersPrec.xml parametersSolver.xml
    EXEDEPS navierStokes_2D_3D_initialSolution)
TRIBITS_COPY_FILES_TO_BINARY_DIR(parameters_navierStokes_2D_3D_initialSolution
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES parametersProblem.xml
    EXEDEPS navierStokes_2D_3D_initialSolution)
TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_navierStokes_2D_3D_initialSolution
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES BFS2d_1600.mesh BFS3dCC.mesh
    EXEDEPS navierStokes_2D_3D_initialSolution)

FOREACH(dimension 2 3)
    FOREACH(initial_mpi_procs 4 6)
        IF(NOT TPL_ENABLE_MPI OR MPI_EXEC_MAX_NUMPROCS LESS initial_mpi_procs OR MPI_EXEC_MAX_NUMPROCS LESS 4)
            CONTINUE()
        ENDIF()
        SET(initial_test_name ${PACKAGE_NAME}_navierStokes_2D_3D_initialSolution_${dimension}D_4_to_${initial_mpi_procs})
        IF(dimension EQUAL 2)
            SET(initial_test_mesh BFS2d_1600.mesh)
        ELSE()
            SET(initial_test_mesh BFS3dCC.mesh)
        ENDIF()
        CONFIGURE_FILE(runInitialSolutionTest.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/${initial_test_name}.cmake @ONLY)
        ADD_TEST(NAME ${initial_test_name}
            COMMAND ${CMAKE_COMMAND} "-DEXECUTABLE=$<TARGET_FILE:${PACKAGE_NAME}_navierStokes_2D_3D_initialSolution>"
                -P ${CMAKE_CURRENT_BINARY_DIR}/${initial_test_name}.cmake)
        SET_TESTS_PROPERTIES(${initial_test_name} PROPERTIES PROCESSORS ${initial_mpi_procs} TIMEOUT 600
            LABELS ${PACKAGE_NAME} REQUIRED_FILES "$<TARGET_FILE:${PACKAGE_NAME}_navierStokes_2D_3D_initialSolution>")
    ENDFOREACH()
ENDFOREACH()
