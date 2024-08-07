TRIBITS_ADD_EXECUTABLE(
	laplaceUnitTest
	SOURCES laplaceUnitTest.cpp
)

TRIBITS_ADD_TEST(
	laplaceUnitTest
	NAME laplace_2D_P1
	ARGS "--dim=2 --FEType=P1"
	NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
	laplaceUnitTest
	NAME laplace_2D_P2
	ARGS "--dim=2 --FEType=P2"
	NUM_MPI_PROCS 4
)
TRIBITS_ADD_TEST(
	laplaceUnitTest
	NAME laplace_3D_P1
	ARGS "--dim=3 --FEType=P1"
	NUM_MPI_PROCS 8
)
TRIBITS_ADD_TEST(
	laplaceUnitTest
	NAME laplace_3D_P2
	ARGS "--dim=3 --FEType=P2"
	NUM_MPI_PROCS 8
)


TRIBITS_COPY_FILES_TO_BINARY_DIR(data_laplaceUnitTest
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrec_Laplace.xml parametersSolver.xml
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrec_Laplace.xml  parametersSolver.xml
	EXEDEPS laplaceUnitTest
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_tests
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ReferenceSolutions
	SOURCE_FILES solution_laplace_2d_P1_4cores.h5 solution_laplace_2d_P2_4cores.h5 solution_laplace_3d_P1_8cores.h5 solution_laplace_3d_P2_8cores.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}/ReferenceSolutions
	DEST_FILES solution_laplace_2d_P1_4cores.h5 solution_laplace_2d_P2_4cores.h5 solution_laplace_3d_P1_8cores.h5 solution_laplace_3d_P2_8cores.h5
	EXEDEPS laplaceUnitTest
)