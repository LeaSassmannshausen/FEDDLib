TRIBITS_ADD_EXECUTABLE(
    fsi_restart
    SOURCES main.cpp
)

# Restart test, two runs of the same case: the uninterrupted run
# (parametersProblemFSI.xml) writes checkpoints at t = 0.01 and 0.02 to
# restartCheckpoints/; the restart (parametersProblemFSI_restart.xml) continues
# from the checkpoint at t = 0.01 and must reproduce the uninterrupted solution
# at t = 0.02.
TRIBITS_ADD_ADVANCED_TEST(
    fsi_restart
    TEST_0 EXEC fsi_restart
        ARGS --problemfile=parametersProblemFSI.xml
        NUM_MPI_PROCS 4
    TEST_1 EXEC fsi_restart
        ARGS --problemfile=parametersProblemFSI_restart.xml
        NUM_MPI_PROCS 4
    OVERALL_NUM_MPI_PROCS 4
    FAIL_FAST
)

SET(fsi_restart_parameter_files
    parametersProblemFSI.xml
    parametersProblemFSI_restart.xml
    parametersPrecFluidMono.xml
    parametersPrecStructure.xml
    parametersPrecGeometry.xml
    parametersSolverFSI.xml
    parametersSolverGeometry.xml
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fsi_restartTest
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${fsi_restart_parameter_files}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${fsi_restart_parameter_files}
    EXEDEPS fsi_restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_fsi_restartTest
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES dfg_fsi_fluid_h004.mesh dfg_fsi_solid_h004.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES dfg_fsi_fluid_h004.mesh dfg_fsi_solid_h004.mesh
    EXEDEPS fsi_restart
)
