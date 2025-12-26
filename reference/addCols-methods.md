# Add columns to LP problem

Add columns (a.k.a. variables/reactions) to an
[LPproblem](LPproblem-class.md).

## Usage

``` r
addCols(lp, ...)

# S4 method for class 'LPproblem_glpk'
addCols(lp, ncols)
```

## Arguments

- lp:

  Object of class [LPproblem](LPproblem-class.md)

- ...:

  Additional parameters passed on to the specific method instance.

- ncols:

  Number of columns to add
