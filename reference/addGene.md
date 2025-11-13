# Add genes or update their data

The function allows you to add one or more genes to a model. When
providing the ID of an already existing genes, you can use this function
to update the gene's information.

## Usage

``` r
addGene(
  model,
  id,
  name = NA,
  CVTerms = NA,
  SBOTerm = rep("SBO:0000243", length(id))
)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- id:

  Character vector with gene IDs

- name:

  Character vector for gene names

- CVTerms:

  Character vector for the genes' CV-Terms

- SBOTerm:

  Character vector for the genes' SBO-Terms
