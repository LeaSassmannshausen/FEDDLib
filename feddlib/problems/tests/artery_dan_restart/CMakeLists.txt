TRIBITS_ADD_EXECUTABLE(
	artery_dan_restart
    SOURCES main.cpp
)

TRIBITS_ADD_TEST(
    artery_dan_restart
    NAME artery_dan_restart
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_artery_dan_restart
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} 
	SOURCE_FILES materialParameters_dan.xml preconditionerParameters_Chemistry.xml 
	preconditionerParameters_Structure.xml simulationParameters.xml solverParameters.xml slurm-script.sh
	referenceSolution/History0.500000.h5
	referenceSolution/History4.000000.h5
	referenceSolution/History25.000000.h5
	referenceSolution/Solution0.h5
	referenceSolution/Solution1.h5
	# referenceSolution/c.h5
	# referenceSolution/c.xmf
	# referenceSolution/d_s.h5
	# referenceSolution/d_s.xmf
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES materialParameters_dan.xml preconditionerParameters_Chemistry.xml 
	preconditionerParameters_Structure.xml simulationParameters.xml solverParameters.xml slurm-script.sh
	History0.500000.h5
	History4.000000.h5
	History25.000000.h5
	Solution0.h5
	Solution1.h5
	# referenceSolution/c.h5
	# referenceSolution/c.xmf
	# referenceSolution/d_s.h5
	# referenceSolution/d_s.xmf
	EXEDEPS artery_dan_restart 
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_artery_dan_restart
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
    SOURCE_FILES Artery_dan_SCI.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES Artery_dan_SCI.mesh
	EXEDEPS artery_dan_restart
)
