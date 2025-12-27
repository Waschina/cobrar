# Flux Balance Analysis

Performs basic flux balance analysis (fba)

## Usage

``` r
fba(model)
```

## Arguments

- model:

  Model of class
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

## Value

A list with flux predictions (reaction fluxes 'fluxes', reduced costs
'redCosts'), and optimization status ()

## See also

Other Flux prediction algorithms:
[`fva()`](https://waschina.github.io/cobrar/reference/fva.md),
[`pfba()`](https://waschina.github.io/cobrar/reference/pfba.md),
[`pfbaHeuristic()`](https://waschina.github.io/cobrar/reference/pfbaHeuristic.md)

## Examples

``` r
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)

# aerobic growth
res_aero <- fba(mod)
cat(" Growth rate:       ", res_aero@obj,"\n",
    "Acetate production:", res_aero@fluxes[mod@react_id == "EX_ac_e"],"\n")
#>  Growth rate:        0.8739215 
#>  Acetate production: 0 

mod <- changeBounds(mod, react = "EX_o2_e", lb = 0) # before: -1000
res_anaero <- fba(mod)
cat(" Growth rate:       ", res_anaero@obj,"\n",
    "Acetate production:", res_anaero@fluxes[mod@react_id == "EX_ac_e"],"\n")
#>  Growth rate:        0.2116629 
#>  Acetate production: 8.503585 
```
