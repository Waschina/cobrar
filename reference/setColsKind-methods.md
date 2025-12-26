# Set column types

Set column/variable types of an [LPproblem](LPproblem-class.md).

## Usage

``` r
setColsKind(lp, ...)

# S4 method for class 'LPproblem_glpk'
setColsKind(lp, j, kind)
```

## Arguments

- lp:

  Object of class [LPproblem](LPproblem-class.md)

- ...:

  Additional parameters passed on to the specific method instance.

- j:

  Indices of columns

- kind:

  Type of respective columns
