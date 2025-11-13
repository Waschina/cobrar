# Couple reaction flux bounds to biomass production

Resets the absolute flux bounds of reactions proportional to the flux
through the biomass reaction.

## Usage

``` r
fluxBMCoupling(
  model,
  BMreact = guessBMReaction(model),
  cpl_c = 400,
  cpl_u = 0.01
)
```

## Arguments

- model:

  Community model of class [ModelComm](ModelComm-class.md)

- BMreact:

  IDs of individual biomass reactions

- cpl_c:

  Coupling constraint `c`

- cpl_u:

  Coupling constraint `u`

## Details

The idea of flux coupling in flux balance models of multi-species
communities was first introduced by Heinken et al. (2013). The idea is
to limit the reactions fluxes in order to prevent that a reaction in one
organism is solely used (i.e., carries a non-zero flux) to benefit
another organism in the community and not the organism that produced the
enzyme for the specific reaction. Therefore, new linear constraints are
added to the flux balance model where the absolute reaction flux bounds
of organism *j* is proportional to the biomass formation of organism
*j*. The coupling constraints are defined as:

\$\$-c v\_{b,j} - u \leq v\_{i,j} \leq c v\_{b,j} + u\$\$

where \\v\_{i,j}\\ is the flux through reaction *i* in organism *j*,
\\v\_{b,j}\\ the biomass reaction of organism *j*. *c* and *u* are
parameters for the coupling constraints that define intercept (*u*) and
slope (*c*) (see figure).

![Coupling Constraints](figures/coupling_constraints.svg)

## References

- Heinken A, Sahoo S, Fleming RMT, Thiele I. Systems-level
  characterization of a host-microbe metabolic symbiosis in the
  mammalian gut. Vol. 4, Gut Microbes; 2013.
  <http://dx.doi.org/10.4161/gmic.22370>
