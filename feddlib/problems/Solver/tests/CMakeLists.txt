### Domain
tribits_add_executable(
facsi
SOURCES facsi.cpp
)

tribits_add_test(
facsi
NAME tests_Solver_facsi
ARGS
NUM_MPI_PROCS 2
)

# Runtime tiers from the laptop baseline in doc/TestRuntimes.md.
feddlib_set_test_runtimes(DEFAULT Low)
