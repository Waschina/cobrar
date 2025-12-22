# Add compartments or update their data

The function allows you to add one or more compartments to a model. When
providing the ID of an already existing compartment, you can use this
function to update the compartment's name.

## Usage

``` r
addCompartment(model, id, name = NA)
```

## Arguments

- model:

  Model of class
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

- id:

  Character vector with compartment IDs

- name:

  Character vector for compartment names

## Examples

``` r
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)
mod <- addCompartment(mod, id = "p", name = "periplasm")
```
