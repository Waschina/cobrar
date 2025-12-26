# Check compartment IDs and Indices

Checks whether compartment IDs or indices are part of (valid) for a
specific model.

## Usage

``` r
checkCompartmentId(model, comp)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- comp:

  A character vector specifying the compartment IDs or a integer vector
  providing the compartment indices in the model.

## Value

A logical vector; TRUE if ID/index is valid, FALSE otherwise.
