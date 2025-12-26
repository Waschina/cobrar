# Index of metabolite(s)

Returns the index(es) of specific metabolite(s).

## Usage

``` r
met_pos(model, met)

# S4 method for class 'ModelOrg,character'
met_pos(model, met)

# S4 method for class 'ModelOrg,numeric'
met_pos(model, met)

# S4 method for class 'ModelOrg,missing'
met_pos(model, met)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- met:

  Character vector with metabolite IDs or Integer vector providing
  indexes.

## Details

Returns NA for metabolite IDs not part of the model or if the index is
larger than the number of metabolites in the model.
