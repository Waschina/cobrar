# Check metabolite IDs and Indices

Checks whether metabolite IDs or indices are part of (valid) for a
specific model.

## Usage

``` r
checkMetId(model, met)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- met:

  A character vector specifying the metabolite IDs or a integer vector
  providing the metabolite indices in the model.

## Value

A logical vector; TRUE if ID/index is valid, FALSE otherwise.
