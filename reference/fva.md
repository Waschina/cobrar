# Flux Variability Analysis (FVA)

Perform Flux Variability Analysis with or without relaxed optimality
constraint

## Usage

``` r
fva(model, react = NULL, opt.factor = 1)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- react:

  Character vector of reaction IDs tested for flux variability. If NULL,
  all reactions are tested.

- opt.factor:

  Numeric value \> 0 to define the required fraction of the objective
  function value. E.g. 0.8 sets the constraint, that in the flux
  variability analysis, the objective function value must at least be
  80% of the original optimal value.

## See also

Other Flux prediction algorithms: [`fba()`](fba.md),
[`pfba()`](pfba.md), [`pfbaHeuristic()`](pfbaHeuristic.md)

## Examples

``` r
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)

# Get flux variability for all exchange reactions
fvares <- fva(mod, react = mod@react_id[grepl("^EX_",mod@react_id)],
              opt.factor = 0.9)
fvares
#>          react growth.fraction   min.flux   max.flux
#> 1      EX_ac_e             0.9   0.000000   3.813556
#> 2   EX_acald_e             0.9   0.000000   2.542370
#> 3     EX_akg_e             0.9   0.000000   1.430083
#> 4     EX_co2_e             0.9  15.206526  26.528850
#> 5    EX_etoh_e             0.9   0.000000   2.214323
#> 6     EX_for_e             0.9   0.000000  11.322324
#> 7     EX_fru_e             0.9   0.000000   0.000000
#> 8     EX_fum_e             0.9   0.000000   0.000000
#> 9  EX_glc__D_e             0.9 -10.000000  -9.046611
#> 10 EX_gln__L_e             0.9   0.000000   0.000000
#> 11 EX_glu__L_e             0.9   0.000000   1.271185
#> 12      EX_h_e             0.9  15.777779  27.100103
#> 13    EX_h2o_e             0.9  20.935920  32.258244
#> 14 EX_lac__D_e             0.9   0.000000   2.145125
#> 15 EX_mal__L_e             0.9   0.000000   0.000000
#> 16    EX_nh4_e             0.9  -5.559972  -4.288787
#> 17     EX_o2_e             0.9 -25.619543 -17.992432
#> 18     EX_pi_e             0.9  -3.214895  -2.893406
#> 19    EX_pyr_e             0.9   0.000000   2.542370
#> 20   EX_succ_e             0.9   0.000000   1.674244
```
