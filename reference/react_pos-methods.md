# Index of reaction(s)

Returns the index(es) of specific reaction(s).

## Usage

``` r
react_pos(model, react)

# S4 method for class 'ModelOrg,character'
react_pos(model, react)

# S4 method for class 'ModelOrg,numeric'
react_pos(model, react)

# S4 method for class 'ModelOrg,missing'
react_pos(model, react)
```

## Arguments

- model:

  Model of class
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

- react:

  Character vector with reaction IDs or Integer vector providing
  indexes.

## Details

Returns NA for reaction IDs not part of the model or the index is larger
than the number of reactions in the model.
