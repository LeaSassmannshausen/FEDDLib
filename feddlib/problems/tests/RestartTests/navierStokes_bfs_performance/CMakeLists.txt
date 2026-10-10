TRIBITS_ADD_EXECUTABLE(navierStokes_bfs_performance SOURCES main.cpp)

SET(performance_inputs parametersProblem.xml parametersPrec.xml parametersSolver.xml)
TRIBITS_COPY_FILES_TO_BINARY_DIR(data_navierStokes_bfs_performance
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${performance_inputs}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${performance_inputs}
    EXEDEPS navierStokes_bfs_performance
)

# Short local smoke runs; larger MPI jobs invoke the same executable directly.
IF(TPL_ENABLE_MPI AND MPI_EXEC_MAX_NUMPROCS GREATER_EQUAL 9)
    FOREACH(dimension 2 3)
        SET(performance_test_name ${PACKAGE_NAME}_navierStokes_bfs_performance_${dimension}D_MPI_9)
        SET(performance_script ${CMAKE_CURRENT_BINARY_DIR}/${performance_test_name}.cmake)
        CONFIGURE_FILE(runPerformance.cmake.in ${performance_script} @ONLY)
        ADD_TEST(NAME ${performance_test_name}
            COMMAND ${CMAKE_COMMAND}
                "-DEXECUTABLE=$<TARGET_FILE:${PACKAGE_NAME}_navierStokes_bfs_performance>"
                -P ${performance_script})
        SET_TESTS_PROPERTIES(${performance_test_name} PROPERTIES
            PROCESSORS 9 LABELS "${PACKAGE_NAME};RestartTests;Performance"
            REQUIRED_FILES "$<TARGET_FILE:${PACKAGE_NAME}_navierStokes_bfs_performance>"
            TIMEOUT 600)
    ENDFOREACH()
ENDIF()
