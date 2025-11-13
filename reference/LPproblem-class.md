# Structure of LPproblem Class

A class structure to link LP problem C++ object.

## Slots

- `ptr`:

  External pointer to LP problem C++ object

- `solver`:

  Solver used for the LP problem

- `method`:

  Specific algorithm used by the LP solver

- `tol_bnd`:

  Numeric value determining how closely the solution must satisfy the
  bounds on variables.
