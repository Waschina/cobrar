# Remove reactions from a model

This function removes specified reactions from a model.

## Usage

``` r
rmReact(model, react, rm_met = TRUE)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- react:

  A character vector stating the reaction IDs in a model or a numeric
  vector providing the reaction indices.

- rm_met:

  Logical. Should metabolites, which are singletons after the reaction
  removal, be deleted as well?

## Value

An updated model of class [ModelOrg](ModelOrg-class.md)

## Note

If the reaction participates in a user constraint, this constraint is
removed from the model.
