# Print reaction(s)

Print the equations of reactions.

## Usage

``` r
printReaction(model, react, use.ids = FALSE)
```

## Arguments

- model:

  Model of class
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

- react:

  A character vector specifying the reaction IDs or a integer vector
  providing the reaction indices in the model.

- use.ids:

  Boolean. Indicating whether metabolite IDs should be printed instead
  of metabolite names.

## Value

A character vector with the individual reaction equations.

## Examples

``` r
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)
# print reaction specified by index
printReaction(mod, react = 2)
#> [1] "(1) Pyruvate + (1) Coenzyme A --> (1) Acetyl-CoA + (1) Formate"
# print reaction specified by ID
printReaction(mod, react = "PFL")
#> [1] "(1) Pyruvate + (1) Coenzyme A --> (1) Acetyl-CoA + (1) Formate"
# print reaction with metabolite IDs instead of metabolite names
printReaction(mod, react = "PFL", use.ids = TRUE)
#> [1] "(1) pyr_c + (1) coa_c --> (1) accoa_c + (1) for_c"
# print multiple reactions at once
printReaction(mod, react = c(2,8))
#> [1] "(1) Pyruvate + (1) Coenzyme A --> (1) Acetyl-CoA + (1) Formate"
#> [2] "(1) D-Glycerate 2-phosphate <==> (1) 3-Phospho-D-glycerate"    
```
