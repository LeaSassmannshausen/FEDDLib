TRIBITS_ADD_EXECUTABLE(
    unsteadyNavierStokes_Restart
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
        TRIBITS_ADD_ADVANCED_TEST(
            unsteadyNavierStokes_Restart_${dimension}D_${mpi_procs}
            OVERALL_WORKING_DIRECTORY TEST_NAME
            TEST_0 COPY_FILES_TO_TEST_DIR
                ${navier_stokes_restart_parameter_files}
                ${navier_stokes_restart_mesh_files}
                SOURCE_DIR ${CMAKE_CURRENT_BINARY_DIR}
            TEST_1 EXEC unsteadyNavierStokes_Restart
                ARGS --problemfile=parametersProblem${dimension}D.xml
                NUM_MPI_PROCS ${mpi_procs}
            TEST_2 EXEC unsteadyNavierStokes_Restart
                ARGS --problemfile=parametersProblem${dimension}D.xml
                     --restartfile=parametersProblem_restart.xml
                NUM_MPI_PROCS ${mpi_procs}
            OVERALL_NUM_MPI_PROCS ${mpi_procs}
            FAIL_FAST
        )
    ENDFOREACH()
ENDFOREACH()

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_unsteadyNavierStokes_Restart
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${navier_stokes_restart_parameter_files}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${navier_stokes_restart_parameter_files}
    EXEDEPS unsteadyNavierStokes_Restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_unsteadyNavierStokes_Restart
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES ${navier_stokes_restart_mesh_files}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${navier_stokes_restart_mesh_files}
    EXEDEPS unsteadyNavierStokes_Restart
)
