# Add metabolites or update their data

The function allows you to add one or more metabolites to a model. When
providing the ID of an already existing metabolite, you can use this
function to update metabolite information.

## Usage

``` r
addMetabolite(
  model,
  id,
  comp = NA,
  name = NA,
  chemicalFormula = NA,
  charge = NA,
  CVTerms = NA,
  SBOTerm = rep("SBO:0000247", length(id))
)
```

## Arguments

- model:

  Model of class
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

- id:

  Character vector with metabolite IDs

- comp:

  Character vector of the metabolites' compartment IDs

- name:

  Character vector for metabolite names

- chemicalFormula:

  Character vector for the metabolites' chemical formulas

- charge:

  Numeric vector for the metabolites' charge

- CVTerms:

  Character vector for the metabolites' CV-Terms

- SBOTerm:

  Character vector for the metabolites' SBO-Terms
