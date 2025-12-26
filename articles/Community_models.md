# Community models

``` r
library(cobrar)
#> Loading required package: Matrix
#> cobrar uses...
#>  - libSBML (v. 5.20.2)
#>  - glpk (v. 5.0)
```

## Merge models of organisms to construct community metabolic models

### Background

Microbial cells in multi-species communities frequently engage in
metabolite exchange interactions[¹](#fn1).

### Workflow

First, we load the wild type *E. coli* core metabolic model.

``` r
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)
```

Next, two single-gene knockout genotypes are generated. For the example
here, we chose the two genes *b0720* (citrate synthase, cs) and *b0767*
(6-phosphogluconolactonase, pgl).

``` r
mod_dCS  <- rmGene(mod, "b0720")
mod_dPGL <- rmGene(mod, "b0767")
```

Now we are ready to generate a community metabolic model by joining the
two models `mod_dCS` and `mod_dPGL`.

``` r
cmod <- joinModels(list(mod_dCS, mod_dPGL))
cmod
#> model ID:                   cobrar-model 
#> model name:                 cobrar-model 
#> number of compartments:     5 
#>                             M1_e  ( M1 extracellular space )
#>                             M1_c  ( M1 cytosol )
#>                             M2_e  ( M2 extracellular space )
#>                             M2_c  ( M2 cytosol )
#>                             e  ( shared extracellular space )
#> number of reactions:        208 
#> number of metabolites:      164 
#> number of unique genes:     272 
#> number of user constraints: 0 
#> number of subsystems:       0 
#> 
#> objective function:         +1 M1_BIOMASS_Ecoli_core_w_GAM +1 M2_BIOMASS_Ecoli_core_w_GAM 
#> objective direction:        maximize
```

The merged model contains five compartments: an extracellular space and
a cytosol for each organism model M1 (“Δcs”) and M2 (“Δpgl”)), plus a
new shared extracellular space `e` via which both genotypes can in
principle exchange metabolites (see figure).

![Figure: Commumnity model](community.svg)

**Figure**: Commumnity model

Before we predict reaction and exchange fluxes for the community, we add
some additional constraints to the community model, namely flux coupling
constraints. These constraints have the effect, that the absolute
reaction flux bounds of organism x are proportional to the biomass
formation of organism x. For more information on flux coupling, please
see the documentation of the function
[`fluxBMCoupling()`](../reference/fluxBMCoupling.md) and the paper by
Heinken et al. (2013)[²](#fn2).

``` r
cmod <- fluxBMCoupling(cmod, cpl_c = 200, cpl_u = 0.01)
```

All set. We can now predict the flux distribution within the community
metabolic model, e.g. via the heuristic parsimonious FBA and store the
results in the variable `sol`.

``` r
sol <- pfbaHeuristic(cmod)
```

To inspect the predicted growth rate per genotype, we can run:

``` r
BMrxn <- guessBMReaction(cmod) # identify the biomass reactions within the community
growth <- sol@fluxes[react_pos(cmod,BMrxn)]
names(growth) <- BMrxn
growth
#> M1_BIOMASS_Ecoli_core_w_GAM M2_BIOMASS_Ecoli_core_w_GAM 
#>                   0.4533128                   0.3778828
```

Please note, that the organisms within the community are numbered
consecutively in the same order as the models are supplied to the
`joinModels(...)` function.

### Analysis of metabolic interchange

Finally, we are interested in the metabolic interactions between both
genotypes. To extract cross-feeding information from the flux
prediction, we can inspect the flux rates of the organism exchange
reactions (see figure above). Metabolites, which have non-zero flux
rates for both respective organism exchange reactions and where the
signs of predicted fluxes are different, are metabolites that are
exchanged between both community members.

First, a data frame is constructed that contains all organism exchange
reaction IDs, names, predicted fluxes, and the corresponding organism
Index (M1 or M2).

``` r
dfCF <- data.frame(exrxn = cmod@react_id[grepl("^M[1|2]_EX_", cmod@react_id)])
dfCF$model <- substr(dfCF$exrxn,1,2)
dfCF$rxnname <- cmod@react_name[react_pos(cmod, dfCF$exrxn)]
dfCF$flux <- sol@fluxes[react_pos(cmod, dfCF$exrxn)]
```

Next we split the data frame by organism:

``` r
dfCF_M1 <- dfCF[dfCF$model == "M1", c("rxnname","flux")]
dfCF_M2 <- dfCF[dfCF$model == "M2", c("rxnname","flux")]
```

Finally, both data frames are merged in a way, that each row represents
a metabolite and with two flux columns, one for the exchange reaction
rates for each genotype. The re-organised data frame makes it easy to
spot cross-fed metabolites.

``` r
dfCF <- merge(dfCF_M1, dfCF_M2, by = "rxnname", suffixes = c(".dCS",".dPGL"))

dfCF
#>                    rxnname   flux.dCS     flux.dPGL
#> 1  2-Oxoglutarate exchange  -1.299013  1.299013e+00
#> 2    Acetaldehyde exchange   0.000000  0.000000e+00
#> 3         Acetate exchange   5.224706 -5.224706e+00
#> 4         Ammonia exchange  -2.471824 -2.060519e+00
#> 5             CO2 exchange  11.864378  1.276368e+01
#> 6      D-Fructose exchange   0.000000  0.000000e+00
#> 7       D-Glucose exchange  -5.851613 -4.148387e+00
#> 8       D-lactate exchange   0.000000  0.000000e+00
#> 9         Ethanol exchange   0.000000  0.000000e+00
#> 10        Formate exchange   0.000000  0.000000e+00
#> 11       Fumarate exchange   0.000000  0.000000e+00
#> 12             H+ exchange  11.720134  4.953649e+00
#> 13            H2O exchange  12.568463  1.811436e+01
#> 14    L-Glutamate exchange   0.000000  8.881784e-16
#> 15    L-Glutamine exchange   0.000000  0.000000e+00
#> 16       L-Malate exchange   0.000000  0.000000e+00
#> 17             O2 exchange -10.041290 -1.362582e+01
#> 18      Phosphate exchange  -1.667602 -1.390117e+00
#> 19       Pyruvate exchange   0.000000  0.000000e+00
#> 20      Succinate exchange   0.000000  0.000000e+00
```

The output shows, that the “Δcs” genotype produces acetate and consumes
2-oxoglutarate, while the “Δpgl” consumes acetate and produced
2-oxoglutarate.

------------------------------------------------------------------------

1.  D’Souza, G., Shitut, S., Preussger, D., Yousif, G., Waschina, S. &
    Kost, C. Ecology and evolution of metabolic cross-feeding
    interactions in bacteria. Nat. Prod. Rep. 35, 455–488 (2018).
    <https://doi.org/10.1039/C8NP00009C>

2.  Heinken, A., Sahoo, S., Fleming, R. M. T. & Thiele, I. Systems-level
    characterization of a host-microbe metabolic symbiosis in the
    mammalian gut. Gut Microbes 4, 28–40 (2013).
    <http://dx.doi.org/10.4161/gmic.22370>
