# Add single constraint

Add a single constraint to an existing [LPproblem](LPproblem-class.md).

## Usage

``` r
addSingleConstraint(lp, ...)

# S4 method for class 'LPproblem_glpk'
addSingleConstraint(lp, coeffs, lb, ub, type)
```

## Arguments

- lp:

  Object of class [LPproblem](LPproblem-class.md)

- ...:

  Additional parameters passed on to the specific method instance.

- coeffs:

  Linear coefficients for variables

- lb:

  Lower bound of constraint

- ub:

  Upper bound of constraint

- type:

  Constraint type
