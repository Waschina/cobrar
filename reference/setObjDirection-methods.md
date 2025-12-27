# Set objective direction

Set the objective function direction in an
[LPproblem](https://waschina.github.io/cobrar/reference/LPproblem-class.md).

## Usage

``` r
setObjDirection(lp, ...)

# S4 method for class 'LPproblem_glpk'
setObjDirection(lp, lpdir)
```

## Arguments

- lp:

  Object of class
  [LPproblem](https://waschina.github.io/cobrar/reference/LPproblem-class.md)

- ...:

  Additional parameters passed on to the specific method instance.

- lpdir:

  Objective direction ("max" or "min")
