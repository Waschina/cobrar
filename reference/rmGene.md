# Remove genes from a model

This function removes specified genes from a model, and optionally also
reactions and metabolites inaccessible after gene knock outs.

## Usage

``` r
rmGene(model, gene, rm_react = TRUE, rm_met = TRUE)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- gene:

  A character vector stating the reaction IDs in a model or a numeric
  vector providing the reaction indices.

- rm_react:

  Logical. Should reaction, which are inaccessible after the gene knock
  outs, be deleted as well?

- rm_met:

  Logical. Should metabolites, which are singletons after the reaction
  removal, be deleted as well?

## Value

An updated model of class [ModelOrg](ModelOrg-class.md)

## Examples

``` r
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)
mod
#> model ID:                   e_coli_core 
#> model name:                 e_coli_core 
#> number of compartments:     2 
#>                             e  ( extracellular space )
#>                             c  ( cytosol )
#> number of reactions:        95 
#> number of metabolites:      72 
#> number of unique genes:     137 
#> number of user constraints: 0 
#> number of subsystems:       0 
#> 
#> objective function:         +1 BIOMASS_Ecoli_core_w_GAM 
#> objective direction:        maximize 

# create a double gene knock-out mutant
mod_KO <- rmGene(mod, c("b4152","b0116"))
mod_KO
#> model ID:                   e_coli_core 
#> model name:                 e_coli_core 
#> number of compartments:     2 
#>                             e  ( extracellular space )
#>                             c  ( cytosol )
#> number of reactions:        92 
#> number of metabolites:      72 
#> number of unique genes:     135 
#> number of user constraints: 0 
#> number of subsystems:       0 
#> 
#> objective function:         +1 BIOMASS_Ecoli_core_w_GAM 
#> objective direction:        maximize 
```
