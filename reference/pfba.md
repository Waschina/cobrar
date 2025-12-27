# Parsimonious Flux Balance Analysis (pFBA)

Performs parsimonious FBA as describe by Lewis et al. 2010.

## Usage

``` r
pfba(model, costcoeffw = NULL, costcoefbw = NULL)
```

## Arguments

- model:

  Model of class
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

- costcoeffw, costcoefbw:

  A numeric vector containing cost coefficients for all
  variables/reactions (forward direction: 'costcoeffw'; backward
  direction: 'costcoefbw'). If set to NULL, all cost coefficients are
  set to 1, so that all variables have the same impact on the objective
  function.

## Value

Returned reduced costs ('redCosts') correspond to the optimization of
the initial linear program (LP), which is basically the initial FBA to
calculate the optimal value of the objective function that is defined in
'model'.

## References

N. E. Lewis et al., “Omic data from evolved E. coli are consistent with
computed optimal growth from genome‐scale models,” Molecular Systems
Biology, vol. 6, no. 1. EMBO, Jan. 2010. doi: 10.1038/msb.2010.47.

## See also

[`pfbaHeuristic()`](https://waschina.github.io/cobrar/reference/pfbaHeuristic.md)

Other Flux prediction algorithms:
[`fba()`](https://waschina.github.io/cobrar/reference/fba.md),
[`fva()`](https://waschina.github.io/cobrar/reference/fva.md),
[`pfbaHeuristic()`](https://waschina.github.io/cobrar/reference/pfbaHeuristic.md)
