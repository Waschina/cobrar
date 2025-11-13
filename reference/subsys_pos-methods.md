# Index of subsystem(s)

Returns the index(es) of specific subsystem(s).

## Usage

``` r
subsys_pos(model, subsys)

# S4 method for class 'ModelOrg,character'
subsys_pos(model, subsys)

# S4 method for class 'ModelOrg,numeric'
subsys_pos(model, subsys)

# S4 method for class 'ModelOrg,missing'
subsys_pos(model, subsys)

# S4 method for class 'ModelOrg,logical'
subsys_pos(model, subsys)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- subsys:

  Character vector with subsystem IDs or Integer vector providing
  indexes.

## Details

Returns NA for subsystem IDs not part of the model or if the index is
larger than the number of subsystems in the model.
