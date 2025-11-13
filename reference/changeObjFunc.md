# Change the objective function

Changes the objective function of a model.

## Usage

``` r
changeObjFunc(model, react, obj_coef = rep(1, length(react)))
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- react:

  Character vector containing the model's reactions IDs, that are part
  of the new objective function

- obj_coef:

  Numeric vector with the objective coefficients of the reactions in the
  vector 'react'.

## Value

An updated model of class [ModelOrg](ModelOrg-class.md)

## Details

Reactions not listed in 'react' are assigned an objective coefficient of
zero.
