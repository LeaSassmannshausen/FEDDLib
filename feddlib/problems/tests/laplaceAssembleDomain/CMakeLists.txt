TRIBITS_ADD_EXECUTABLE(
	laplaceAssembleDomain
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	laplaceAssembleDomain
    NAME laplaceAssembleDomain_DIM_2_P1
    ARGS "--DIM=2"
    NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
	laplaceAssembleDomain
    NAME laplaceAssembleDomain_DIM_3_P1
    ARGS "--DIM=3"
    NUM_MPI_PROCS 1
)
TRIBITS_ADD_TEST(
	laplaceAssembleDomain
    NAME laplaceAssembleDomain_DIM_3_P1
    ARGS "--DIM=3"
    NUM_MPI_PROCS 8
)
TRIBITS_ADD_TEST(
	laplaceAssembleDomain
    NAME laplaceAssembleDomain_DIM_2_P2
    ARGS "--DIM=2 --FEType=P2"
    NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
	laplaceAssembleDomain
    NAME laplaceAssembleDomain_DIM_3_P2
    ARGS "--DIM=3 --FEType=P2"
    NUM_MPI_PROCS 1
)
TRIBITS_ADD_TEST(
	laplaceAssembleDomain
    NAME laplaceAssembleDomain_DIM_3_P2
    ARGS "--DIM=3 --FEType=P2"
    NUM_MPI_PROCS 8
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_laplaceAssembleDomain
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml 
	EXEDEPS laplaceAssembleDomain
)
