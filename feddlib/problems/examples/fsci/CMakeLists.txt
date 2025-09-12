TRIBITS_ADD_EXECUTABLE(
	fsci
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fsci
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml parametersPrecChem.xml parametersProblemFSCI.xml parametersPrecStructureCE.xml parametersSolverFSCI.xml parametersSolverGeometry.xml parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecStructure.xml  laplace_parabolic_parabolic_fsi_fluid_2mm_P2.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml parametersPrecChem.xml  parametersProblemFSCI.xml parametersPrecStructureCE.xml parametersSolverFSCI.xml parametersSolverGeometry.xml parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecStructure.xml  laplace_parabolic_parabolic_fsi_fluid_2mm_P2.h5
)
TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_fsci
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES fsi_fluid_2mm.mesh fsi_solid_2mm.mesh
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES fsi_fluid_2mm.mesh fsi_solid_2mm.mesh
	EXEDEPS fsci
)	