# Count elements in formulas

Counts the number of elements in molecules based on their chemical
formulas.

## Usage

``` r
countElements(formulas, multiplier = 1)
```

## Arguments

- formulas:

  Character vector of chemical formulas.

- multiplier:

  Multiplier for the element coutns in each molecule

## Value

A numeric matrix with named columns for individual elements. Rows
contain the counts of elements for each formula in the same order as
formulas were provided to the function.

## Details

Formulas may not have element numbers containing decimal points.

## Examples

``` r
countElements(c("C6H12O6","C48H72CoN11O8","HCOOH"))
#>       C  H O Co  N
#> [1,]  6 12 6  0  0
#> [2,] 48 72 8  1 11
#> [3,]  1  2 2  0  0
```
