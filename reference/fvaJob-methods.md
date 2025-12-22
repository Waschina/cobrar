# Wrapper function for efficient FVA

Wrapper function for efficient FVA

## Usage

``` r
fvaJob(lp, ...)

# S4 method for class 'LPproblem_glpk'
fvaJob(lp, ind)
```

## Arguments

- lp:

  Object of class
  [LPproblem](https://waschina.github.io/cobrar/reference/LPproblem-class.md)

- ...:

  Additional parameters passed on to the specific method instance.

- ind:

  Indices of variables to be tested in FVA
