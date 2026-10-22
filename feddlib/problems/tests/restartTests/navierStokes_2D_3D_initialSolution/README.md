# Navier–Stokes initial solution tests

The 2D `BFS2d_1600.mesh` and 3D `BFS3dCC.mesh` cases initialize on four or six
MPI ranks from a converged steady solution computed on four ranks.

Each isolated test runs these phases:

1. Solve the steady Navier–Stokes problem and check nonlinear convergence. Save
   only velocity, pressure and compatibility metadata at source time `0.015`.
2. Initialize a reference by directly assigning the saved vectors. Evolve it
   from time zero to `0.008` with the new viscosity and integration settings.
3. Load the initial solution through the production initialization path. Check
   exact field equality, time zero and empty history, including after solver
   setup. A zero-duration run must leave the fields and history unchanged.
4. Run the initialized simulation and compare the written states at `0`, `0.004`
   and `0.008` against direct initialization. Print absolute and relative l2 errors
   for velocity and pressure; require exact equality at zero and relative error
   at most `1e-12` thereafter.

The source metadata uses BDF1 and `dt=0.0025`; the new simulation uses BDF2 and
`dt=0.004`. The first new timestep uses the normal BDF1 startup. No source history
exists, and the source time is off the new timestep grid. The 2D four-rank case
also checks rejection of a mesh with a modified coordinate.

The mesh overrides and linear solver settings are shared with
`../navierStokes_2D_3D_bfs`. Initial solution settings are in `parametersProblem.xml`.

Build and run from the workspace root:

```sh
cmake --build build --target problems_navierStokes_2D_3D_initialSolution
ctest --test-dir build --output-on-failure -R navierStokes_2D_3D_initialSolution
```
