# Reuse the existing 2D BFS2d_1600 Navier–Stokes case and restart executable.
if(TPL_ENABLE_MPI AND MPI_EXEC_MAX_NUMPROCS GREATER_EQUAL 4)
    set(recovery_model NavierStokes)
    set(recovery_executable ${PACKAGE_NAME}_navierStokes_2D_3D_bfs)
    set(recovery_case_directory ${CMAKE_CURRENT_SOURCE_DIR}/../navierStokes_2D_3D_bfs)
    set(recovery_input_files parametersProblem.xml parametersProblem2D.xml parametersPrec.xml parametersSolver.xml)
    set(recovery_mesh_files BFS2d_1600.mesh)
    set(recovery_test_name ${PACKAGE_NAME}_navierStokes_2D_3D_bfs_recovery)
    configure_file(../runRecoveryTest.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/${recovery_test_name}.cmake @ONLY)
    add_test(NAME ${recovery_test_name}
        COMMAND ${CMAKE_COMMAND} "-DEXECUTABLE=$<TARGET_FILE:${recovery_executable}>"
            -P ${CMAKE_CURRENT_BINARY_DIR}/${recovery_test_name}.cmake)
    set_tests_properties(${recovery_test_name} PROPERTIES PROCESSORS 4 TIMEOUT 600 LABELS ${PACKAGE_NAME}
        REQUIRED_FILES "$<TARGET_FILE:${recovery_executable}>")
endif()
