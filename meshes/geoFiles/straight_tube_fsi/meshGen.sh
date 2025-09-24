gmsh -3 -optimize geometry.geo
gmsh -3 fluid.geo
gmsh -3 solid.geo
#./gmsh2medit.py fluid.msh fluid.mesh
#./gmsh2medit.py solid.msh solid.mesh
