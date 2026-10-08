TRIBITS_ADD_EXECUTABLE(unsteadyNonLinElasticity_restart SOURCES main.cpp)

SET(elasticity_restart_inputs
    parametersProblem.xml parametersProblem_restart.xml parametersPrec.xml parametersSolver.xml)

SET(reference_mpi_procs 4)
FOREACH(restart_mpi_procs 4 6)
    IF(NOT TPL_ENABLE_MPI OR reference_mpi_procs GREATER MPI_EXEC_MAX_NUMPROCS
        OR restart_mpi_procs GREATER MPI_EXEC_MAX_NUMPROCS)
        CONTINUE()
    ENDIF()
    SET(restart_test_name ${PACKAGE_NAME}_unsteadyNonLinElasticity_restart_2D_${reference_mpi_procs}_to_${restart_mpi_procs})
    SET(restart_test_script ${CMAKE_CURRENT_BINARY_DIR}/${restart_test_name}.cmake)
    CONFIGURE_FILE(runRestartTest.cmake.in ${restart_test_script} @ONLY)
    ADD_TEST(NAME ${restart_test_name}
        COMMAND ${CMAKE_COMMAND}
            "-DEXECUTABLE=$<TARGET_FILE:${PACKAGE_NAME}_unsteadyNonLinElasticity_restart>"
            -P ${restart_test_script})
    SET_TESTS_PROPERTIES(${restart_test_name} PROPERTIES
        PROCESSORS ${restart_mpi_procs}
        LABELS ${PACKAGE_NAME}
        REQUIRED_FILES "$<TARGET_FILE:${PACKAGE_NAME}_unsteadyNonLinElasticity_restart>")
ENDFOREACH()

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadyNonLinElasticity_restart
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${elasticity_restart_inputs}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${elasticity_restart_inputs}
    EXEDEPS unsteadyNonLinElasticity_restart)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadyNonLinElasticity_restart
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES square_h02.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES square_h02.mesh
    EXEDEPS unsteadyNonLinElasticity_restart)
