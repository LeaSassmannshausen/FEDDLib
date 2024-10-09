TRIBITS_ADD_EXECUTABLE(
	fsci_restart
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fsci_restart
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml parametersPrecChem.xml parametersProblemFSCI.xml parametersProblemFluid.xml parametersProblemStructure.xml parametersSolverFSCI.xml parametersSolverGeometry.xml parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecStructure.xml parametersPrecStructureCE.xml plistProblemLaplace.xml plistPrecLaplace.xml plistSolverLaplace.xml 
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml parametersPrecChem.xml  parametersProblemFSCI.xml parametersProblemFluid.xml parametersProblemStructure.xml parametersPrecStructureCE.xml parametersSolverFSCI.xml parametersSolverGeometry.xml parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecStructure.xml plistProblemLaplace.xml plistPrecLaplace.xml plistSolverLaplace.xml p
	EXEDEPS fsci_restart
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_fsci_restart
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES fluid_length_0_5_simp.mesh solid_length_0_5_simp.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES  fluid_length_0_5_simp.mesh solid_length_0_5_simp.mesh
	EXEDEPS fscfsci_restarti
)
