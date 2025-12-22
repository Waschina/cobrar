# Identify dead end metabolites

Searches a metabolic network for metabolites that can be produced but
not consumed, and vice versa.

## Usage

``` r
deadEndMetabolites(object)
```

## Arguments

- object:

  Model of class
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

## Value

A list with two elements: "dem" is a character vector with the IDs of
identified dead end metabolites. "der" is a character vector with IDs of
reactions, that involve at least one of the dead end metabolites.

## Note

The algorithm is adapted from the original 'sybil' package, which in
turn is adapted from the implementation in 'cobratoolbox'. Important
detail: The algorithm function in the 'sybil' package may return
different results than this cobrar implementation. This is because the
'sybil' function does not take into account the possibility that an
irreversible reaction is defined in the direction RHS to LHS (lower
bound \> 0 and upper bound = 0).

## Examples

``` r
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)
deadEndMetabolites(mod)
#> $dem
#> [1] "gln__L_e" "mal__L_e" "fru_e"    "fum_e"   
#> 
#> $der
#> [1] "EX_fru_e"    "EX_fum_e"    "EX_gln__L_e" "EX_mal__L_e" "FRUpts2"    
#> [6] "FUMt2_2"     "GLNabc"      "MALt2_2"    
#> 
```
