TRIBITS_ADD_EXECUTABLE(
	sci_Cube_Restart
    SOURCES main.cpp
)

TRIBITS_ADD_TEST(
    sci_Cube_Restart
    NAME sci_Cube_Restart
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_sci_Cube_Restart
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecChem.xml parametersPrecCE.xml parametersProblemSCI.xml parametersProblemStructure.xml parametersSolverSCI.xml parametersPrecStructure.xml slurm-script.sh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecChem.xml parametersPrecCE.xml parametersProblemSCI.xml parametersProblemStructure.xml parametersSolverSCI.xml parametersPrecStructure.xml slurm-script.sh
	EXEDEPS sci_Cube_Restart
)

# TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_sciTest
# 	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes/SPP2311
#     SOURCE_FILES artery_h_025_L_025.mesh artery_h_0_8_L_0_25.mesh plaque_short_6_k.mesh
# 	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
# 	DEST_FILES artery_h_025_L_025.mesh artery_h_0_8_L_0_25.mesh plaque_short_6_k.mesh
# 	EXEDEPS sci_Cube_Restart
# )
