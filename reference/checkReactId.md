# Check reaction IDs and Indices

Checks whether reaction IDs or indices are part of (valid) for a
specific model.

## Usage

``` r
checkReactId(model, react)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- react:

  A character vector specifying the reaction IDs or a integer vector
  providing the reaction indices in the model.

## Value

A logical vector; TRUE if ID/index is valid, FALSE otherwise.
