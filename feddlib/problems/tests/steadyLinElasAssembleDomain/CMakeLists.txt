TRIBITS_ADD_EXECUTABLE(
	steadyLinElasAssembleDomain
	SOURCES main.cpp
)


TRIBITS_ADD_TEST(
	steadyLinElasAssembleDomain
    NAME steadyLinElasAssembleDomain_DIM_2_P1
    ARGS "--DIM=2 --FEType=P1 --MeshName=square.mesh"
    NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
	steadyLinElasAssembleDomain
    NAME steadyLinElasAssembleDomain_DIM_2_P2
    ARGS "--DIM=2 --FEType=P2 --MeshName=square.mesh"
    NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
	steadyLinElasAssembleDomain
    NAME steadyLinElasAssembleDomain_DIM_3_P1
    ARGS "--DIM=3 --FEType=P1 --MeshName=cube_h_1.mesh"
    NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
	steadyLinElasAssembleDomain
    NAME steadyLinElasAssembleDomain_DIM_3_P2
    ARGS "--DIM=3 --FEType=P2 --MeshName=cube_h_1.mesh"
    NUM_MPI_PROCS 4
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(data_steadyLinElasAssembleDomain
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	EXEDEPS steadyLinElasAssembleDomain
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(mesh_steadyLinElasAssembleDomain
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
	SOURCE_FILES cube_h_1.mesh square.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES cube_h_1.mesh square.mesh
	EXEDEPS steadyLinElasAssembleDomain
)
