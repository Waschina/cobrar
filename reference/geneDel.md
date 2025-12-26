# Identify reactions affected by gene knockouts

This function identifies reactions, which cannot be catalyzed anymore
when a specified set of genes is deleted from a model.

## Usage

``` r
geneDel(model, gene)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- gene:

  Character or numeric vector providing the IDs or indices of genes to
  be deleted from 'model'.

## Value

Character vector with reactions IDs.
