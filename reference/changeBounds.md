# Change flux bounds

The function changes either upper bounds, lower bounds, or both for
specific reactions.

## Usage

``` r
changeBounds(model, react, lb = NULL, ub = NULL)
```

## Arguments

- model:

  Model of class
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)

- react:

  A character vector stating the reaction IDs in a model or a numeric
  vector providing the reaction indices.

- lb:

  A numeric vector giving the new lower flux bounds for reactions
  `react`. If `lb` is of length 1, the same value will be used for all
  reactions.

- ub:

  A numeric vector giving the new upper flux bounds for reactions
  `react`. If `ub` is of length 1, the same value will be used for all
  reactions.

## Value

An updated model of class
[ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)
