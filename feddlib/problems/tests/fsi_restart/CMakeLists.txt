TRIBITS_ADD_EXECUTABLE(
	fsi_restart
    SOURCES main.cpp
)

TRIBITS_ADD_TEST(
    fsi_restart
    NAME fsi_restart
    ARGS
    NUM_MPI_PROCS 2
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fsi_restartTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml parametersProblemFSI.xml parametersProblemFluid.xml parametersProblemStructure.xml parametersSolverFSI.xml parametersSolverGeometry.xml parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecStructure.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml  parametersProblemFSI.xml parametersProblemFluid.xml parametersProblemStructure.xml parametersSolverFSI.xml parametersSolverGeometry.xml parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecStructure.xml
	EXEDEPS fsi_restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_fsi_restartTest
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES dfg_fsi_fluid_h004.mesh dfg_fsi_solid_h004.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES dfg_fsi_fluid_h004.mesh dfg_fsi_solid_h004.mesh
	EXEDEPS fsi_restart
)
