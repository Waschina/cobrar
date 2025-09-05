#' Structure of LPproblem Class
#'
#' A class structure to link LP problem C++ object.
#'
#' @slot ptr External pointer to LP problem C++ object
#' @slot solver Solver used for the LP problem
#' @slot method Specific algorithm used by the LP solver
#' @slot tol_bnd Numeric value determining how closely the solution must
#' satisfy the bounds on variables.
#'
#' @aliases LPproblem
#'
#' @exportClass LPproblem
setClass("LPproblem",
         slots = c(
           ptr = "externalptr",
           solver = "character",
           method = "character",
           tol_bnd = "numeric"
         )
)


# Set generic functions, where each solver gets its own version

#' Initialize a LP problem
#'
#' Transfers variables, constraints and objectives to \link{LPproblem}.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname loadLPprob-methods
#' @export
setGeneric("loadLPprob", function(lp, ...) {
  standardGeneric("loadLPprob")
})

#' Set objective direction
#'
#' Set the objective function direction in an \link{LPproblem}.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname setObjDirection-methods
#' @export
setGeneric("setObjDirection", function(lp, ...) {
  standardGeneric("setObjDirection")
})

#' Add columns to LP problem
#'
#' Add columns (a.k.a. variables/reactions) to an \link{LPproblem}.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname addCols-methods
#' @export
setGeneric("addCols", function(lp, ...) {
  standardGeneric("addCols")
})

#' Add rows to LP problem
#'
#' Add rows (a.k.a. constraints) to an \link{LPproblem}.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname addRows-methods
#' @export
setGeneric("addRows", function(lp, ...) {
  standardGeneric("addRows")
})

#' Populate a constraint-X-variable matrix
#'
#' Add linear coefficients to the constraint-X-variable matrix of an
#' \link{LPproblem}.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname loadMatrix-methods
#' @export
setGeneric("loadMatrix", function(lp, ...) {
  standardGeneric("loadMatrix")
})

#' Set column bounds and objective coefficients
#'
#' Set column bounds and objective coefficients of an \link{LPproblem}.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname setColsBndsObjCoefs-methods
#' @export
setGeneric("setColsBndsObjCoefs", function(lp, ...) {
  standardGeneric("setColsBndsObjCoefs")
})

#' Set column types
#'
#' Set column/variable types of an \link{LPproblem}.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname setColsKind-methods
#' @export
setGeneric("setColsKind", function(lp, ...) {
  standardGeneric("setColsKind")
})

#' Set row bounds
#'
#' Set row bounds of an \link{LPproblem}.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname setRowsBnds-methods
#' @export
setGeneric("setRowsBnds", function(lp, ...) {
  standardGeneric("setRowsBnds")
})

#' Solve an LP problem
#'
#' Solves an \link{LPproblem}.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname solveLp-methods
#' @export
setGeneric("solveLp", function(lp, ...) {
  standardGeneric("solveLp")
})

#' Get the objective value of a solved LP problem
#'
#' Get the objective value of a solved \link{LPproblem}.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname getObjValue-methods
#' @export
setGeneric("getObjValue", function(lp, ...) {
  standardGeneric("getObjValue")
})

#' Get the solver status
#'
#' Get the solver status of a solved \link{LPproblem}.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname getSolStat-methods
#' @export
setGeneric("getSolStat", function(lp, ...) {
  standardGeneric("getSolStat")
})

#' Retrieve column primal value
#'
#' Retrieve column primal value of a solved \link{LPproblem}.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname getColsPrimal-methods
#' @export
setGeneric("getColsPrimal", function(lp, ...) {
  standardGeneric("getColsPrimal")
})

#' Retrieve column reduced costs
#'
#' Retrieve column reduced costs of a solved \link{LPproblem}.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname getRedCosts-methods
#' @export
setGeneric("getRedCosts", function(lp, ...) {
  standardGeneric("getRedCosts")
})

#' Add single constraint
#'
#' Add a single constraint to an existing \link{LPproblem}.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname addSingleConstraint-methods
#' @export
setGeneric("addSingleConstraint", function(lp, ...) {
  standardGeneric("addSingleConstraint")
})

#' Wrapper function for efficient FVA
#'
#' Wrapper function for efficient FVA
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname fvaJob-methods
#' @export
setGeneric("fvaJob", function(lp, ...) {
  standardGeneric("fvaJob")
})


#' Delete an LP problem
#'
#' Deletes an existing \link{LPproblem} and frees associated memory.
#'
#' @param lp Object of class \link{LPproblem}
#' @param ... Additional parameters passed on to the specific method instance.
#'
#' @docType methods
#' @rdname deleteLP-methods
#' @export
setGeneric("deleteLP", function(lp, ...) {
  standardGeneric("deleteLP")
})

