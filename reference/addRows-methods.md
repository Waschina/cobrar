# Add rows to LP problem

Add rows (a.k.a. constraints) to an
[LPproblem](https://waschina.github.io/cobrar/reference/LPproblem-class.md).

## Usage

``` r
addRows(lp, ...)

# S4 method for class 'LPproblem_glpk'
addRows(lp, nrows)
```

## Arguments

- lp:

  Object of class
  [LPproblem](https://waschina.github.io/cobrar/reference/LPproblem-class.md)

- ...:

  Additional parameters passed on to the specific method instance.

- nrows:

  Number of rows to add
