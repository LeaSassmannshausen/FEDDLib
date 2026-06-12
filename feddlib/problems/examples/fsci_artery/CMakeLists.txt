TRIBITS_ADD_EXECUTABLE(
	fsci_artery
    SOURCES main.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(data_fsciArtery
	SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
	SOURCE_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml 
	parametersProblemFSCI.xml parametersSolverFSCI.xml parametersSolverGeometry.xml 
	parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecFluidBlock.xml 
	parametersPrecStructure.xml
	 laplace_parabolic_arterie_with_plaque_initma_media_adventitia_fluid_fine_P2.h5
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES parametersPrecGE.xml parametersPrecGI.xml parametersPrecGeometry.xml  
	parametersProblemFSCI.xml parametersSolverFSCI.xml parametersSolverGeometry.xml 
	parametersPrecFluidMono.xml parametersPrecFluidTeko.xml parametersPrecFluidBlock.xml 
	parametersPrecStructure.xml 
	laplace_parabolic_arterie_with_plaque_initma_media_adventitia_fluid_fine_P2.h5
	EXEDEPS fsci_artery
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(meshes_fsciArtery
	SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES arterie_with_plaque_initma_media_adventitia_fluid_fine.mesh arterie_with_plaque_initma_media_adventitia_solid_fine.mesh  
	DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
	DEST_FILES arterie_with_plaque_initma_media_adventitia_fluid_fine.mesh arterie_with_plaque_initma_media_adventitia_solid_fine.mesh 
	EXEDEPS fsci_artery
)
