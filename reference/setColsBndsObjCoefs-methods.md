# Set column bounds and objective coefficients

Set column bounds and objective coefficients of an
[LPproblem](https://waschina.github.io/cobrar/reference/LPproblem-class.md).

## Usage

``` r
setColsBndsObjCoefs(lp, ...)

# S4 method for class 'LPproblem_glpk'
setColsBndsObjCoefs(lp, j, lb, ub, obj_coef, type = NULL)
```

## Arguments

- lp:

  Object of class
  [LPproblem](https://waschina.github.io/cobrar/reference/LPproblem-class.md)

- ...:

  Additional parameters passed on to the specific method instance.

- j:

  Indices of variables to be updated

- lb:

  New lower bound of variable

- ub:

  New upper bound of variable

- obj_coef:

  New object coefficient of variable

- type:

  New type of column/variable
