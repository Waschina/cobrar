# Check subsystem IDs and Indices

Checks whether subsystem IDs or indices are part of (valid) for a
specific model.

## Usage

``` r
checkSubsystemId(model, subsystem)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- subsystem:

  A character vector specifying the subsystem IDs or a integer vector
  providing the subsystem indices in the model.

## Value

A logical vector; TRUE if ID/index is valid, FALSE otherwise.
