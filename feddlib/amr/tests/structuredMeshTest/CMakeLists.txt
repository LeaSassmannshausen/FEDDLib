TRIBITS_ADD_EXECUTABLE(
	structuredMeshTest
	SOURCES main.cpp
)

TRIBITS_ADD_TEST(
	structuredMeshTest
    NAME structuredMeshTest2D
    ARGS 
	--problemfile=parametersProblem.xml
    NUM_MPI_PROCS 4
)

TRIBITS_ADD_TEST(
	structuredMeshTest
    NAME structuredMeshTest3D
    ARGS 
	--problemfile=parametersProblem3D.xml
    NUM_MPI_PROCS 4
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_structuredMeshTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml parametersProblem3D.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec.xml  parametersProblem.xml parametersSolver.xml parametersProblem3D.xml 
	EXEDEPS structuredMeshTest
)


