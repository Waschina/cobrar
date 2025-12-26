# Index of compartment(s)

Returns the index(es) of specific compartment(s).

## Usage

``` r
comp_pos(model, comp)

# S4 method for class 'ModelOrg,character'
comp_pos(model, comp)

# S4 method for class 'ModelOrg,numeric'
comp_pos(model, comp)

# S4 method for class 'ModelOrg,missing'
comp_pos(model, comp)

# S4 method for class 'ModelOrg,logical'
comp_pos(model, comp)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- comp:

  Character vector with compartment IDs or Integer vector providing
  indexes.

## Details

Returns NA for compartment IDs not part of the model or if the index is
larger than the number of compartments in the model.
