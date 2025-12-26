# Introduction to cobrar

``` r
library(cobrar)
#> Loading required package: Matrix
#> cobrar uses...
#>  - libSBML (v. 5.20.2)
#>  - glpk (v. 5.0)
```

## Performing Flux Balance Analysis (FBA)

The following (very brief) example illustrates how to use `cobrar` to
perform flux balance analysis for a core metabolism model of
*Escherichia coli*. In a first simulation, fluxes and growth are
predicted for the default constraints of the model, which represents a
minimal medium with glucose as the sole carbon source and the presence
of oxygen. The second simulation performs the same flux balance analysis
but without oxygen to simulate an anoxic growth environment.

First, we load the wild type *E. coli* core metabolic model.

``` r
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)
```

Next, FBA is performed to predict the model’s growth rate and reaction
fluxes.

``` r
sol <- fba(mod)
sol
#> Algorithm:              FBA 
#> Solver status:          solution is optimal 
#> Optimization status:    optimization process was successful 
#> Objective fct. value:   0.8739215 
#> Secondary objective:    NA
```

We can also inspect the metabolite consumption and production by looking
at the predicted fluxes for all exchange reactions.

``` r
getExchanges(mod, sol)
#>             ID                    name       flux
#> 1      EX_ac_e        Acetate exchange   0.000000
#> 2   EX_acald_e   Acetaldehyde exchange   0.000000
#> 3     EX_akg_e 2-Oxoglutarate exchange   0.000000
#> 4     EX_co2_e            CO2 exchange  22.809833
#> 5    EX_etoh_e        Ethanol exchange   0.000000
#> 6     EX_for_e        Formate exchange   0.000000
#> 7     EX_fru_e     D-Fructose exchange   0.000000
#> 8     EX_fum_e       Fumarate exchange   0.000000
#> 9  EX_glc__D_e      D-Glucose exchange -10.000000
#> 10 EX_gln__L_e    L-Glutamine exchange   0.000000
#> 11 EX_glu__L_e    L-Glutamate exchange   0.000000
#> 12      EX_h_e             H+ exchange  17.530865
#> 13    EX_h2o_e            H2O exchange  29.175827
#> 14 EX_lac__D_e      D-lactate exchange   0.000000
#> 15 EX_mal__L_e       L-Malate exchange   0.000000
#> 16    EX_nh4_e        Ammonia exchange  -4.765319
#> 17     EX_o2_e             O2 exchange -21.799493
#> 18     EX_pi_e      Phosphate exchange  -3.214895
#> 19    EX_pyr_e       Pyruvate exchange   0.000000
#> 20   EX_succ_e      Succinate exchange   0.000000
```

This table indicates that the organism is growing aerobically, because
oxygen is consumed (reaction “EX_o2_e” has a flux below 0).

To simulate anaerobic growth, the lower bound of the exchange reaction
“EX_o2_e” can be set to 0.

``` r
mod_anaero <- changeBounds(mod, react = "EX_o2_e", lb = 0)
```

Now, the growth rate and reaction fluxes can be predicted for anaerobic
growth.

``` r
sol_anaero <- fba(mod_anaero)
sol_anaero
#> Algorithm:              FBA 
#> Solver status:          solution is optimal 
#> Optimization status:    optimization process was successful 
#> Objective fct. value:   0.2116629 
#> Secondary objective:    NA
getExchanges(mod_anaero, sol_anaero)
#>             ID                    name          flux
#> 1      EX_ac_e        Acetate exchange  8.503585e+00
#> 2   EX_acald_e   Acetaldehyde exchange  0.000000e+00
#> 3     EX_akg_e 2-Oxoglutarate exchange  0.000000e+00
#> 4     EX_co2_e            CO2 exchange -3.781782e-01
#> 5    EX_etoh_e        Ethanol exchange  8.279455e+00
#> 6     EX_for_e        Formate exchange  1.780467e+01
#> 7     EX_fru_e     D-Fructose exchange  0.000000e+00
#> 8     EX_fum_e       Fumarate exchange  0.000000e+00
#> 9  EX_glc__D_e      D-Glucose exchange -1.000000e+01
#> 10 EX_gln__L_e    L-Glutamine exchange  0.000000e+00
#> 11 EX_glu__L_e    L-Glutamate exchange  0.000000e+00
#> 12      EX_h_e             H+ exchange  3.055422e+01
#> 13    EX_h2o_e            H2O exchange -7.115796e+00
#> 14 EX_lac__D_e      D-lactate exchange  0.000000e+00
#> 15 EX_mal__L_e       L-Malate exchange  0.000000e+00
#> 16    EX_nh4_e        Ammonia exchange -1.154156e+00
#> 17     EX_o2_e             O2 exchange  0.000000e+00
#> 18     EX_pi_e      Phosphate exchange -7.786445e-01
#> 19    EX_pyr_e       Pyruvate exchange  0.000000e+00
#> 20   EX_succ_e      Succinate exchange -9.215424e-14
```

## Editing a metabolic network model

This example adds the 4-aminobutyrate degradation pathway to the *E.
coli* core metabolic model, which was already loaded above and is stored
in the object named `mod`.

First, the transporter reaction is added to the model, that can import
4-aminobutyrate (ID: `4abut`).

``` r
mod <- addReact(mod, id = "ABUTt", Scoef = c(-1,-1,1,1),
                met = c("4abut_e","h_e","4abut_c","h_c"), reversible = TRUE,
                lb = -1000, ub = 1000,
                reactName = "4-aminobutyrate transport in via proton symport",
                metName = c("4-aminobutyrate",NA, "4-aminobutyrate",NA),
                metComp = c("e","e","c","c"), metCharge = c(0,NA,0,NA),
                metChemicalFormula = c("C4H9NO2",NA,"C4H9NO2",NA))
```

Next, we add the exchange reaction for 4-aminobutyrate.

``` r
mod <- addReact(mod, id = "EX_4abut_e", Scoef = c(-1), met = "4abut_e",
                lb = -1.5, ub = 1000, reactName = "4-aminobutyrate exchange")
```

The lower bound of this exchange reaction is set to -1.5 mmol/(gDW\*hr)
to simulate the availability of the metabolite.

As a next step, the reaction 4-aminobutyrate amninotransferase (EC
[2.6.1.19](https://www.rhea-db.org/rhea/23352)) is added to the model,
which transfers the amino group of 4-aminobutyrate to
alpha-ketoglutarate (ID `akg_c`) and forms L-glutamate (ID `glu__L_c`)
and succinic semialdehyde (ID `sucsal_c`).

``` r
mod <- addReact(mod, id = "ABTA", Scoef = c(-1,-1,1,1),
                met = c("4abut_c","akg_c","glu__L_c","sucsal_c"),
                lb = 0,
                reactName = "4-aminobutyrate transaminase",
                metName = c(NA,NA,NA,"Succinic semialdehyde"),
                metComp = c(NA,NA,NA,"c"), metCharge = c(NA,NA,NA,-1),
                metChemicalFormula = c(NA,NA,NA,"C4H5O3"),
                subsystem = "GABAdegr", subsystemName = "4-aminobutyrate degradation",
                CVTerms = "bqbiol_is;http://identifiers.org/ec-code/2.6.1.19",
                gprAssoc = "b2662 | b1302")
```

Finally, the reaction Succinate-semialdehyde dehydrogenase (EC 1.2.1.24)
is added, which oxidises succinic semialdehyde to form succinate (ID
`succ_c`).

``` r
mod <- addReact(mod, id = "SSALx", Scoef = c(-1,-1,-1,2,1,1),
                met = c("h2o_c","nad_c","sucsal_c","h_c","nadh_c","succ_c"),
                lb = 0,
                reactName = "Succinate-semialdehyde dehydrogenase (NAD)",
                subsystem = "GABAdegr",
                CVTerms = "bqbiol_is;http://identifiers.org/ec-code/1.2.1.24",
                gprAssoc = "b1525")
```

All new reaction can be printed as equations:

``` r
printReaction(mod, c("EX_4abut_e","ABUTt","ABTA","SSALx"))
#> [1] "(1) 4-aminobutyrate <==> "                                                               
#> [2] "(1) H+ + (1) 4-aminobutyrate <==> (1) H+ + (1) 4-aminobutyrate"                          
#> [3] "(1) 2-Oxoglutarate + (1) 4-aminobutyrate --> (1) L-Glutamate + (1) Succinic semialdehyde"
#> [4] "(1) H2O + (1) NAD + (1) Succinic semialdehyde --> (2) H+ + (1) NADH + (1) Succinate"
```

By performing an FBA on the updated model, we can see that the growth
rate increased with the new pathway and that 4-aminobutyrate is indeed
consumed by *E. coli*.

``` r
sol <- fba(mod)
sol
#> Algorithm:              FBA 
#> Solver status:          solution is optimal 
#> Optimization status:    optimization process was successful 
#> Objective fct. value:   0.9674959 
#> Secondary objective:    NA
getExchanges(mod, sol)
#>             ID                     name       flux
#> 1      EX_ac_e         Acetate exchange   0.000000
#> 2   EX_acald_e    Acetaldehyde exchange   0.000000
#> 3     EX_akg_e  2-Oxoglutarate exchange   0.000000
#> 4     EX_co2_e             CO2 exchange  24.827727
#> 5    EX_etoh_e         Ethanol exchange   0.000000
#> 6     EX_for_e         Formate exchange   0.000000
#> 7     EX_fru_e      D-Fructose exchange   0.000000
#> 8     EX_fum_e        Fumarate exchange   0.000000
#> 9  EX_glc__D_e       D-Glucose exchange -10.000000
#> 10 EX_gln__L_e     L-Glutamine exchange   0.000000
#> 11 EX_glu__L_e     L-Glutamate exchange   0.000000
#> 12      EX_h_e              H+ exchange  17.907968
#> 13    EX_h2o_e             H2O exchange  30.375354
#> 14 EX_lac__D_e       D-lactate exchange   0.000000
#> 15 EX_mal__L_e        L-Malate exchange   0.000000
#> 16    EX_nh4_e         Ammonia exchange  -3.775562
#> 17     EX_o2_e              O2 exchange -24.459205
#> 18     EX_pi_e       Phosphate exchange  -3.559127
#> 19    EX_pyr_e        Pyruvate exchange   0.000000
#> 20   EX_succ_e       Succinate exchange   0.000000
#> 21  EX_4abut_e 4-aminobutyrate exchange  -1.500000
```
