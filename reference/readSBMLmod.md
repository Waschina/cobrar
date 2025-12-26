# Reads an SBML file and constructs an object of class 'ModelOrg'

Reads an SBML file and constructs an object of class 'ModelOrg'

## Usage

``` r
readSBMLmod(file_path)
```

## Arguments

- file_path:

  Path to SBML file.

## Value

A [ModelOrg-class](ModelOrg-class.md) object.

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
```
