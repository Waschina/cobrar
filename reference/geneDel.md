# Identify reactions affected by gene knockouts

This function identifies reactions, which cannot be catalyzed anymore
when a specified set of genes is deleted from a model.

## Usage

``` r
geneDel(model, gene, single = FALSE)
```

## Arguments

- model:

  Model of class
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

- gene:

  Character or numeric vector providing the IDs or indices of genes to
  be deleted from 'model'.

- single:

  If FALSE (default), the effect of simultaneous gene deletions of all
  genes in `gene` are assumed. If TRUE, results for single gene
  deletions are returned.

## Value

Character vector with reactions IDs if `single` is FALSE, an a list of
character vectors if `single` is TRUE. In the latter case, the i-th
element in the returned list corresponds to the gene deletion results of
the i-th gene of the input parameter.

## Examples

``` r
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)

# Identify reactions that would be lost after multiple gene deletions
geneDel(mod, gene = c("b3916","b1723","b4025"))
#> [1] "PFK" "PGI"

# Identify reactions that would be lost after single gene deletions
geneDel(mod, gene = c("b3916","b1723","b4025"), single = TRUE)
#> [[1]]
#> character(0)
#> 
#> [[2]]
#> character(0)
#> 
#> [[3]]
#> [1] "PGI"
#> 
```
