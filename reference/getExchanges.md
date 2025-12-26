# Get metabolite exchange rates

Summarize the predicted fluxes for exchange reactions.

## Usage

``` r
getExchanges(model, sol)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- sol:

  Flux prediction results as an object of class
  [FluxPrediction](FluxPrediction-class.md).

## Value

A data.frame that summarizes the predicted metabolite exchange rates
(=fluxes of exchange reactions).
