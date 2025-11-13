# Add constraints to model

Add linear reaction flux constraints to a metabolic network.

## Usage

``` r
addConstraint(model, react, coeff, rtype, lb = NULL, ub = NULL)

# S4 method for class 'ModelOrg,character,numeric,character'
addConstraint(model, react, coeff, rtype, lb = NULL, ub = NULL)

# S4 method for class 'ModelOrg,list,list,character'
addConstraint(model, react, coeff, rtype, lb = NULL, ub = NULL)
```

## Arguments

- model:

  Model of class [ModelOrg](ModelOrg-class.md)

- react:

  Character vector or a list of character vectors containing the model's
  reactions IDs that are part of the respective constraint.

- coeff:

  Numeric vector or list of numeric vectors defining the coefficients
  for the reactions listed in 'react'.

- rtype:

  Character vector describing the type of the constraint(s). See
  details.

- lb, ub:

  Numeric vector defining the lower and upper bound(s) of the
  constraint(s).

## Details

The slot "rtype" describes the type of each constraint. Valid values and
their effects are:

|        |                                    |                            |
|--------|------------------------------------|----------------------------|
| *code* | *description*                      | *rule*                     |
| "F"    | free constraint                    | \\-\infty \< x \< \infty\\ |
| "L"    | constraint with lower bound        | \\lb \leq x \leq \infty\\  |
| "U"    | constraint with upper bound        | \\-\infty \leq x \leq ub\\ |
| "D"    | double-bounded (ranged) constraint | \\lb \leq x \leq ub\\      |
| "E"    | fixed (equality constraint)        | \\lb = x = ub\\            |

## Examples

``` r
fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
mod <- readSBMLmod(fpath)

# Simulate anaerobic growth
mod <- changeBounds(mod, "EX_o2_e", lb = 0)

# Limit the proton production depending on the growth rate
mod <- addConstraint(mod,
                     react = c("EX_h_e","BIOMASS_Ecoli_core_w_GAM"),
                     coeff = c(1, -20), rtype = "U", ub = 0)

```
