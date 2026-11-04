tribits_add_executable(
    fsi_artery_boundary_tests
    SOURCES main.cpp
)

set(boundary_test_parameters
    parametersPrecGE.xml
    parametersPrecGI.xml
    parametersPrecGeometry.xml
    parametersProblemFSI.xml
    parametersSolverFSI.xml
    parametersSolverGeometry.xml
    parametersPrecFluidMono.xml
    parametersPrecFluidTeko.xml
    parametersPrecFluidBlock.xml
    parametersPrecStructure.xml
)

tribits_copy_files_to_binary_dir(data_fsiArteryBoundaryTests
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${boundary_test_parameters}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${boundary_test_parameters}
    EXEDEPS fsi_artery_boundary_tests
)

# Reuse the original artery inputs until the boundary-test geometry is supplied.
tribits_copy_files_to_binary_dir(profiles_fsiArteryBoundaryTests
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../fsi_artery
    SOURCE_FILES laplace_parabolic_artery_P1.h5 laplace_parabolic_artery_P2.h5
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES laplace_parabolic_artery_P1.h5 laplace_parabolic_artery_P2.h5
    EXEDEPS fsi_artery_boundary_tests
)

tribits_copy_files_to_binary_dir(meshes_fsiArteryBoundaryTests
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES arterie_with_plaque_initma_media_adventitia_fluid_fine.mesh arterie_with_plaque_initma_media_adventitia_solid_fine.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES arterie_with_plaque_initma_media_adventitia_fluid_fine.mesh arterie_with_plaque_initma_media_adventitia_solid_fine.mesh
    EXEDEPS fsi_artery_boundary_tests
)
