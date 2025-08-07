TRIBITS_ADD_EXECUTABLE(
	sci_Test
    SOURCES main.cpp
)

TRIBITS_ADD_TEST(
    sci_Test
    NAME sci_Test
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_sciTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecChem.xml parametersPrecCE.xml parametersProblemSCI.xml parametersSolverSCI.xml parametersPrecStructure.xml slurm-script.sh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecChem.xml parametersPrecCE.xml parametersProblemSCI.xml parametersSolverSCI.xml parametersPrecStructure.xml slurm-script.sh
	EXEDEPS sci_Test
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_sciTest
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/
    SOURCE_FILES fsi_solid_2mm.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES fsi_solid_2mm.mesh
	EXEDEPS sci_Test
)
