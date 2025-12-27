# Fixate biomass ratios

Biomass ratios are fixed via introducing new linear constraints.

## Usage

``` r
fixBMRatios(model, BMreact = guessBMReaction(model), tolerance = 0)
```

## Arguments

- model:

  Community model of class
  [ModelComm](https://waschina.github.io/cobrar/reference/ModelComm-class.md)

- BMreact:

  IDs of individual biomass reactions

- tolerance:

  Tolerated deviation of ratios

## Details

Fixate the biomass ratios in a community models according to models'
relative abundances.
