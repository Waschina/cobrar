# Additional constraints

``` r
library(cobrar)
#> Loading required package: Matrix
#> cobrar uses...
#>  - libSBML (v. 5.20.2)
#>  - glpk (v. 5.0)
```

## Constraint for a multi-carbon source growth environment

Standard FBA models already contain the linear constraints:

- Stationarity constraints: $\mathbf{S}\ \mathbf{v} = \mathbf{0}$
- Flux bounds on individual reactions:
  ${\mathbf{l}\mathbf{b}} \leq \mathbf{v} \leq {\mathbf{u}\mathbf{b}}$

In specific cases, it can be useful to include additional linear
constraints to a model, e.g., to add additional thermodynamic
constraints or interactions between reactions.

This tutorial aims to illustrate how additional constraints can be added
to an existing model and to inspect their impact on flux distribution
predictions.

First, we load a model of *Escherichia coli*’s core metabolism.

``` r
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)
```

Second, we change the lower bounds of the exchange reactions for three
potential carbon sources, namely Glucose, Fructose, and Fumarate. By
doing so, we simulate a growth environment that provides three carbon
sources simultaneously. The maximum uptake rate for each carbon source
is set to $3\ mmol\ gDW^{- 1}\ hr^{- 1}$. Furthermore, we set the lower
bound for the exchange reaction of oxygen (`EX_o2_e`) to
$0\ mmol\ gDW^{- 1}\ hr^{- 1}$, to simulate anoxic conditions.

``` r
mod <- changeBounds(mod, react = c("EX_glc__D_e","EX_fru_e", "EX_fum_e"),
                    lb = c(-3,-3,-3))
mod <- changeBounds(mod, react = "EX_o2_e", lb = 0)
```

Fumarate has 4 carbon atoms, while Glucose and Fructose have 6 carbon
atoms each. Here, we want to add an additional constraint to the model
that limits carbon source uptake to a total of 35 C-atom-mmol per gram
dry weight and and per hour. This can be accomplished by:

``` r
mod <- addConstraint(mod, react = c("EX_glc__D_e","EX_fru_e", "EX_fum_e"),
                     coeff = c(6,6,4), rtype = "L", lb = -35)
# print the user constraint
printConstraint(mod)
#> [1] "-35 <= +6 EX_fru_e +4 EX_fum_e +6 EX_glc__D_e < Inf"
```

Next, we can perform flux balance analysis
([`fba()`](../reference/fba.md)) with the additional constraint and
inspect predicted fluxes of all exchange reactions.

``` r
sol <- fba(mod)
getExchanges(mod, sol)
#>             ID                    name          flux
#> 1      EX_ac_e        Acetate exchange  5.234167e+00
#> 2   EX_acald_e   Acetaldehyde exchange  0.000000e+00
#> 3     EX_akg_e 2-Oxoglutarate exchange  0.000000e+00
#> 4     EX_co2_e            CO2 exchange -1.514230e-01
#> 5    EX_etoh_e        Ethanol exchange  5.144425e+00
#> 6     EX_for_e        Formate exchange  1.078766e+01
#> 7     EX_fru_e     D-Fructose exchange -3.000000e+00
#> 8     EX_fum_e       Fumarate exchange -2.664535e-15
#> 9  EX_glc__D_e      D-Glucose exchange -2.833333e+00
#> 10 EX_gln__L_e    L-Glutamine exchange  0.000000e+00
#> 11 EX_glu__L_e    L-Glutamate exchange  0.000000e+00
#> 12      EX_h_e             H+ exchange  1.772191e+01
#> 13    EX_h2o_e            H2O exchange -4.678495e+00
#> 14 EX_lac__D_e      D-lactate exchange  0.000000e+00
#> 15 EX_mal__L_e       L-Malate exchange  0.000000e+00
#> 16    EX_nh4_e        Ammonia exchange -4.621253e-01
#> 17     EX_o2_e             O2 exchange  0.000000e+00
#> 18     EX_pi_e      Phosphate exchange -3.117702e-01
#> 19    EX_pyr_e       Pyruvate exchange  0.000000e+00
#> 20   EX_succ_e      Succinate exchange  0.000000e+00
```

As wen can see, Fructose and Glucose are consumed, while the uptake of
Fumarate is virtually zero. Also, there is a quite high production of
protons (ID: `EX_h_e`). Next, we will add another constraint, that
limits the release of protons depending on the growth rate.
Specifically, we will limit the release to 5 mmol per newly formed gram
dry weight biomass: $$v_{H^{+}} \leq 5\ v_{Biomass}$$

Which is the same as: $$v_{H^{+}} - 5\ v_{Biomass} \leq 0$$

This constraint can easily be added to the model:

``` r
mod <- addConstraint(mod, react = c("EX_h_e","BIOMASS_Ecoli_core_w_GAM"),
                     coeff = c(1, -5), rtype = "U", ub = 0)
sol <- fba(mod)
getExchanges(mod, sol)
#>             ID                    name       flux
#> 1      EX_ac_e        Acetate exchange  0.0000000
#> 2   EX_acald_e   Acetaldehyde exchange  0.0000000
#> 3     EX_akg_e 2-Oxoglutarate exchange  0.0000000
#> 4     EX_co2_e            CO2 exchange 11.0422250
#> 5    EX_etoh_e        Ethanol exchange 10.4801777
#> 6     EX_for_e        Formate exchange  0.2438012
#> 7     EX_fru_e     D-Fructose exchange -3.0000000
#> 8     EX_fum_e       Fumarate exchange -0.6091407
#> 9  EX_glc__D_e      D-Glucose exchange -2.4272395
#> 10 EX_gln__L_e    L-Glutamine exchange  0.0000000
#> 11 EX_glu__L_e    L-Glutamate exchange  0.0000000
#> 12      EX_h_e             H+ exchange  0.3235326
#> 13    EX_h2o_e            H2O exchange -0.1848859
#> 14 EX_lac__D_e      D-lactate exchange  0.0000000
#> 15 EX_mal__L_e       L-Malate exchange  0.0000000
#> 16    EX_nh4_e        Ammonia exchange -0.3528317
#> 17     EX_o2_e             O2 exchange  0.0000000
#> 18     EX_pi_e      Phosphate exchange -0.2380359
#> 19    EX_pyr_e       Pyruvate exchange  0.0000000
#> 20   EX_succ_e      Succinate exchange  0.0000000
```

The results show that as an effect of the second constraint the proton
production (`EX_h_e`) is reduced, and all three carbon sources are
consumed, including Fumarate.

## Adding multiple constraints at once

In the above example, both constraints were added to the model
successively, each with a separate
[`addConstraint()`](../reference/addConstraint-methods.md) call. You can
achieve the same final outcome model by adding both constraints
simultaneously:

``` r
# Load model
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)

# Change bounds
mod <- changeBounds(mod, react = c("EX_glc__D_e","EX_fru_e", "EX_fum_e"),
                    lb = c(-3,-3,-3))
mod <- changeBounds(mod, react = "EX_o2_e", lb = 0)

# Adding constraints
mod <- addConstraint(mod,
                     react = list(
                       c("EX_glc__D_e","EX_fru_e", "EX_fum_e"),
                       c("EX_h_e","BIOMASS_Ecoli_core_w_GAM")
                     ),
                     coeff = list(
                       c(6,6,4),
                       c(1,-5)
                     ), rtype = c("L","U"), lb = c(-35,NA), ub = c(NA, 0))

# Perform FBA
sol <- fba(mod)
getExchanges(mod, sol)
#>             ID                    name       flux
#> 1      EX_ac_e        Acetate exchange  0.0000000
#> 2   EX_acald_e   Acetaldehyde exchange  0.0000000
#> 3     EX_akg_e 2-Oxoglutarate exchange  0.0000000
#> 4     EX_co2_e            CO2 exchange 11.0422250
#> 5    EX_etoh_e        Ethanol exchange 10.4801777
#> 6     EX_for_e        Formate exchange  0.2438012
#> 7     EX_fru_e     D-Fructose exchange -3.0000000
#> 8     EX_fum_e       Fumarate exchange -0.6091407
#> 9  EX_glc__D_e      D-Glucose exchange -2.4272395
#> 10 EX_gln__L_e    L-Glutamine exchange  0.0000000
#> 11 EX_glu__L_e    L-Glutamate exchange  0.0000000
#> 12      EX_h_e             H+ exchange  0.3235326
#> 13    EX_h2o_e            H2O exchange -0.1848859
#> 14 EX_lac__D_e      D-lactate exchange  0.0000000
#> 15 EX_mal__L_e       L-Malate exchange  0.0000000
#> 16    EX_nh4_e        Ammonia exchange -0.3528317
#> 17     EX_o2_e             O2 exchange  0.0000000
#> 18     EX_pi_e      Phosphate exchange -0.2380359
#> 19    EX_pyr_e       Pyruvate exchange  0.0000000
#> 20   EX_succ_e      Succinate exchange  0.0000000
```
