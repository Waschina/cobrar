# Initialize a LP problem

Transfers variables, constraints and objectives to
[LPproblem](https://waschina.github.io/cobrar/reference/LPproblem-class.md).

## Usage

``` r
loadLPprob(lp, ...)

# S4 method for class 'LPproblem_glpk'
loadLPprob(
  lp,
  nCols,
  nRows,
  mat,
  ub,
  lb,
  obj,
  rlb,
  rtype,
  lpdir,
  rub = NULL,
  ctype = NULL
)
```

## Arguments

- lp:

  Object of class
  [LPproblem](https://waschina.github.io/cobrar/reference/LPproblem-class.md)

- ...:

  Additional parameters passed on to the specific method instance.

- nCols:

  Number of columns/variables

- nRows:

  Number of rows/constraints

- mat:

  constraint-X-variable coefficient matrix

- ub:

  Variable values' upper bounds

- lb:

  Variable values' lower bounds

- obj:

  Linear objective coefficients for variables

- rlb:

  Constraints' lower bounds

- rtype:

  Constraint types

- lpdir:

  Objective direction ("max" or "min")

- rub:

  Constraints' upper bounds

- ctype:

  Variable types
