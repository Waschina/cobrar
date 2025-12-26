# Set row bounds

Set row bounds of an [LPproblem](LPproblem-class.md).

## Usage

``` r
setRowsBnds(lp, ...)

# S4 method for class 'LPproblem_glpk'
setRowsBnds(lp, i, lb, ub, type)
```

## Arguments

- lp:

  Object of class [LPproblem](LPproblem-class.md)

- ...:

  Additional parameters passed on to the specific method instance.

- i:

  Indices of rows/constraints

- lb:

  Lower bounds of constraints

- ub:

  Upper bound of constraints

- type:

  Type of constraint bounds
