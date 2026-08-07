TRIBITS_ADD_EXECUTABLE(
	steadyNonLinElas
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	steadyNonLinElas
	NAME steadyNonLinElas
	ARGS
	NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_steadyNonLinElas
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS steadyNonLinElas
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_steadyNonLinElas
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES fsi_solid_2mm.mesh square_h02.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES fsi_solid_2mm.mesh square_h02.mesh
	EXEDEPS steadyNonLinElas
)
