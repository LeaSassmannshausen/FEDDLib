tribits_add_executable(unsteadyNavierStokes_restart SOURCES main.cpp)

set(restart_example_inputs
    parametersProblem.xml parametersProblem_restart.xml
    parametersPrec.xml parametersSolver.xml README.md)
tribits_copy_files_to_binary_dir(data_unsteadyNavierStokes_restart
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    SOURCE_FILES ${restart_example_inputs}
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES ${restart_example_inputs}
    EXEDEPS unsteadyNavierStokes_restart
)

tribits_copy_files_to_binary_dir(mesh_unsteadyNavierStokes_restart
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/meshes
    SOURCE_FILES BFS3dCC.mesh
    DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
    DEST_FILES BFS3dCC.mesh
    EXEDEPS unsteadyNavierStokes_restart
)
