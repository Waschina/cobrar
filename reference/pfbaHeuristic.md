# Heuristic parsimonious Flux Balance Analysis (pFBA)

Performs a heuristic version of the parsimonious FBA algorithm. See
details.

## Usage

``` r
pfbaHeuristic(model, costcoeffw = NULL, costcoefbw = NULL, pFBAcoeff = 1e-06)
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

- pFBAcoeff:

  Numeric value to weight the minimization of total flux within the
  combined objective function. See details.

## Details

The exact-solution pFBA algorithm described by Lewis et al. 2010
consists of two optimization steps: (1) A basic flux balance analysis is
performed to obtain the optimal value of the objective function (e.g.,
growth rate). (2) The objective function from (1) is transformed into a
constraint where the value of the function is fixed to the determined
optimal value. A new objective is defined that minimizes the absolute
sum of fluxes through the metabolic network.  
The here implemented heuristic pFBA performs only one optimization step.
Therefore, the original objective function is combined with a term that
minimizes the absolute (weighted) sum of fluxes: \$\$max:
\sum\_{i=1}^{n}(c_i v_i) - p \sum\_{i=1}^{n}(w_i \|v_i\|)\$\$ for
maximization direction or \$\$min: \sum\_{i=1}^{n}(c_i v_i) + p
\sum\_{i=1}^{n}(w_i \|v_i\|)\$\$ if the original objective function is
minimized.  
Here, \\c_i\\ is the original objective coefficient of reaction \\i\\,
\\v_i\\ the flux through reaction \\i\\, \\n\\ the number of reactions
in the model, and \\w_i\\ the weight of reaction \\i\\ (see arguments
'costcoeffw', 'costcoefbw'). The scalar parameter \\p\\ (argument
'pFBAcoeff') defines the weighting of the total flux minimization.
Increasing this value increases the weight of total flux minimization,
possibly at costs for the value of the objective function defined in
'model' (e.g., flux through biomass reaction).  
This heuristic implementation of a pFBA is the core of the gap-filling
algorithm of 'gapseq' (Zimmermann et al. 2021).

## References

N. E. Lewis et al., “Omic data from evolved E. coli are consistent with
computed optimal growth from genome‐scale models,” Molecular Systems
Biology, vol. 6, no. 1. EMBO, Jan. 2010. doi: 10.1038/msb.2010.47.

J. Zimmermann, C. Kaleta, and S. Waschina, “gapseq: informed prediction
of bacterial metabolic pathways and reconstruction of accurate metabolic
models,” Genome Biology, vol. 22, no. 1. Springer Science and Business
Media LLC, Mar. 10, 2021. doi: 10.1186/s13059-021-02295-1.

## See also

Other Flux prediction algorithms:
[`fba()`](https://waschina.github.io/cobrar/reference/fba.md),
[`fva()`](https://waschina.github.io/cobrar/reference/fva.md),
[`pfba()`](https://waschina.github.io/cobrar/reference/pfba.md)
