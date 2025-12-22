# Structure of Constraints Class

This class represents user constraints that can be added to a model of
class
[ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)
in addition to the stationarity constraint (\\S v = 0\\) and flux
bounds.

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

## Slots

- `coeff`:

  A sparse numeric matrix of
  [dgCMatrix-class](https://rdrr.io/pkg/Matrix/man/dgCMatrix-class.html)
  representing the coefficients for each reaction in the model. Each row
  denotes a user constraint, each column a reaction in the model in the
  same order as in slot "S" in the corresponding
  [ModelOrg](https://waschina.github.io/cobrar/reference/ModelOrg-class.md)
  object.

- `lb`:

  Numeric vector providing the lower bound for each constraint.

- `ub`:

  Numeric vector providing the lower bound for each constraint.

- `rtype`:

  Character vector stating the constraint type. See details.
