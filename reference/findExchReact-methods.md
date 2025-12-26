# Find exchange reactions

Finds all exchange reactions within a (community) model and report them
in a data frame.

## Usage

``` r
findExchReact(model)

# S4 method for class 'ModelOrg'
findExchReact(model)

# S4 method for class 'ModelComm'
findExchReact(model)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md) or
  [ModelComm](ModelComm-class.md)

## Value

If 'model' is of class [ModelOrg](ModelOrg-class.md), a data.frame is
returned stating all exchange reaction IDs, their index in the reaction
list, the respective metabolite id, name, and index in the metabolite
list. If the model is of class [ModelComm](ModelComm-class.md), the
resulting data.frame contains the community exchange reactions and the
organism-specific exchange reaction. The latter are exchange reactions,
that connect extracellular metabolites of the organism metabolic network
models with the shared extracellular space.
