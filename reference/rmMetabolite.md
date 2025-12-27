# Remove metabolites from a model

This function removes specified metabolites from a model.

## Usage

``` r
rmMetabolite(model, met)
```

## Arguments

- model:

  Model of class
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

- met:

  A character vector stating the metabolite IDs in a model or a numeric
  vector providing the metabolite indices.

## Value

An updated model of class
[ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

## Note

If at least one of the provided metabolites still participates in a
reaction, the function stops with an error message.
