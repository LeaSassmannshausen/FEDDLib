TRIBITS_ADD_EXECUTABLE(
	laplaceAssembleDomain
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	laplaceAssembleDomain
    NAME laplaceAssembleDomain
    ARGS
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_laplaceAssembleDomain
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml 
	EXEDEPS laplaceAssembleDomain
)
