# Populate a constraint-X-variable matrix

Add linear coefficients to the constraint-X-variable matrix of an
[LPproblem](LPproblem-class.md).

## Usage

``` r
loadMatrix(lp, ...)

# S4 method for class 'LPproblem_glpk'
loadMatrix(lp, ne, ia, ja, ra)
```

## Arguments

- lp:

  Object of class [LPproblem](LPproblem-class.md)

- ...:

  Additional parameters passed on to the specific method instance.

- ne:

  Number of nonzero coefficients

- ia:

  row indices of nonzero coefficients

- ja:

  column indices of nonzero coefficients

- ra:

  nonzero values of respective entries
