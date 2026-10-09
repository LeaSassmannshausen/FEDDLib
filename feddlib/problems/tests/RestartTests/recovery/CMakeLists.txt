# Exercise recovery using the existing 2D PDE cases and restart executables.
foreach(recovery_model NavierStokes FSI)
    if(NOT TPL_ENABLE_MPI OR MPI_EXEC_MAX_NUMPROCS LESS 4)
        continue()
    endif()
    if(recovery_model STREQUAL "NavierStokes")
        set(recovery_executable ${PACKAGE_NAME}_unsteadyNavierStokes_restart)
        set(recovery_case_directory ${CMAKE_CURRENT_SOURCE_DIR}/../unsteadyNavierStokes_restart)
        set(recovery_input_files parametersProblem.xml parametersProblem2D.xml parametersPrec.xml parametersSolver.xml)
        set(recovery_mesh_files BFS2d_1600.mesh)
    else()
        set(recovery_executable ${PACKAGE_NAME}_fsi_2D_turek)
        set(recovery_case_directory ${CMAKE_CURRENT_SOURCE_DIR}/../fsi_2D_turek)
        set(recovery_input_files parametersProblemFSI.xml parametersPrecFluidMono.xml parametersPrecStructure.xml
            parametersPrecGeometry.xml parametersSolverFSI.xml parametersSolverGeometry.xml)
        set(recovery_mesh_files dfg_fsi_fluid_h004.mesh dfg_fsi_solid_h004.mesh)
    endif()
    set(recovery_test_name ${PACKAGE_NAME}_recovery_${recovery_model}_2D)
    configure_file(runRecoveryTest.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/${recovery_test_name}.cmake @ONLY)
    add_test(NAME ${recovery_test_name}
        COMMAND ${CMAKE_COMMAND} "-DEXECUTABLE=$<TARGET_FILE:${recovery_executable}>"
            -P ${CMAKE_CURRENT_BINARY_DIR}/${recovery_test_name}.cmake)
    set_tests_properties(${recovery_test_name} PROPERTIES PROCESSORS 4 TIMEOUT 600 LABELS ${PACKAGE_NAME}
        REQUIRED_FILES "$<TARGET_FILE:${recovery_executable}>")
endforeach()
