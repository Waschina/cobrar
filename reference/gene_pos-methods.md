# Index of gene(s)

Returns the index(es) of specific gene(s).

## Usage

``` r
gene_pos(model, gene)

# S4 method for class 'ModelOrg,character'
gene_pos(model, gene)

# S4 method for class 'ModelOrg,numeric'
gene_pos(model, gene)

# S4 method for class 'ModelOrg,missing'
gene_pos(model, gene)
```

## Arguments

- model:

  Model of class
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

- gene:

  Character vector with gene IDs or Integer vector providing indexes.

## Details

Returns NA for gene IDs not part of the model or if the index is larger
than the number of genes in the model.
