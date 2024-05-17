#' Structure of Constraints Class
#'
#' This class represents user constraints that can be added to a model of class
#' \link{ModelOrg} in addition to the stationarity constraint (\eqn{S v = 0})
#' and flux bounds.
#'
#' @slot coeff A sparse numeric matrix of \link[Matrix]{dgCMatrix-class}
#' representing the coefficients for each reaction in the model. Each row
#' denotes a user constraint, each column a reaction in the model in the same
#' order as in slot "S" in the corresponding \link{ModelOrg} object.
#' @slot lb Numeric vector providing the lower bound for each constraint.
#' @slot ub Numeric vector providing the lower bound for each constraint.
#' @slot rtype Character vector stating the constraint type. See details.
#'
#' @details
#' The slot "rtype" describes the type of each constraint. Valid values and
#' their effects are:
#' | *code* | *description* | *rule* |
#' | :----: | :--- | :----: |
#' | "F" | free constraint | \eqn{-\infty < x < \infty} |
#' | "L" | constraint with lower bound | \eqn{lb \leq x \leq \infty} |
#' | "U" | constraint with upper bound | \eqn{-\infty \leq x \leq ub} |
#' | "D" | double-bounded (ranged) constraint | \eqn{lb \leq x \leq ub} |
#' | "E" | fixed (equality constraint) | \eqn{lb = x = ub} |
#'
#' @aliases Constraints
#'
#' @family Object classes
#' @exportClass Constraints
setClass("Constraints",
         slots = c(
           coeff = "dgCMatrix",
           lb = "numeric",
           ub = "numeric",
           rtype = "character"
         )
)
