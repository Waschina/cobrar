# Check gene IDs and Indices

Checks whether gene IDs or indices are part of (valid) for a specific
model.

## Usage

``` r
checkGeneId(model, gene)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- gene:

  A character vector specifying the gene IDs or a integer vector
  providing the gene indices in the model.

## Value

A logical vector; TRUE if ID/index is valid, FALSE otherwise.
