# Get metabolite exchange rates

Summarize the predicted fluxes for exchange reactions.

## Usage

``` r
getExchanges(model, sol)
```

## Arguments

- model:

  Model of class
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

- sol:

  Flux prediction results as an object of class
  [FluxPrediction](https://waschina.github.io/cobrar/reference/FluxPrediction-class.md).

## Value

A data.frame that summarizes the predicted metabolite exchange rates
(=fluxes of exchange reactions).
