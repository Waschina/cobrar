# Add subsystems or update their data

The function allows you to add one or more subsystems to a model. When
providing the ID of an already existing subsystem, you can use this
function to update the subsystem's name.

## Usage

``` r
addSubsystem(model, id, name = NA)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- id:

  Character vector with subsystem IDs

- name:

  Character vector for subsystem names

## Examples

``` r
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)
mod <- addSubsystem(mod, id = "Bifidoshunt",
                    name = "glucose fermentation to acetate and lactate (Bifidobacteria)")
```
