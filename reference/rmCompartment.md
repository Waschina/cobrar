# Remove compartments from a model

This function removes specified compartments from a model.

## Usage

``` r
rmCompartment(model, comp)
```

## Arguments

- model:

  Model of class
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

- comp:

  A character vector stating the compartment IDs in a model or a numeric
  vector providing the compartment indices.

## Value

An updated model of class
[ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

## Note

If at least one of the provided compartments still has metabolites
associated with it, the function stops with an error message.
