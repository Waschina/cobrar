#' Class for flux prediction results
#'
#' This class represents a flux prediction results, e.g., from Flux Balance
#' Analysis or derived methods.
#'
#' @slot algorithm Algorithm used for flux prediction.
#' @slot ok The LP solver's code for the type of optimization result.
#' @slot ok_term  The LP solver's term for the type of optimization result.
#' @slot stat Generic status (integer code) of the current basic solution of the
#' optimization problem.
#' @slot stat_term Generic status (character term) of the current basic solution
#' of the optimization problem.
#' @slot obj Objective value.
#' @slot obj_sec Value of secondary objective function. E.g.: Summed absolute
#' fluxes in \link{pfba}.
#' @slot fluxes Predicted flux values.
#' @slot redCosts Predicted reduced costs (or "*dual value*") for reactions.
#'
#' @aliases FluxPrediction
#'
#' @family Flux analysis tools
#' @exportClass FluxPrediction
setClass("FluxPrediction",
         slots = c(
           algorithm = "character",
           ok = "integer",
           ok_term = "character",
           stat = "integer",
           stat_term = "character",
           obj = "numeric",
           obj_sec = "numeric",
           fluxes = "numeric",
           redCosts = "numeric"
         )
)


#' Print a short summary of a flux prediction result
#'
#' Displays a summary of the results obtained from a metabolic flux prediction
#' (e.g., FBA or pFBA).
#'
#' @param object S4-object of class \link{FluxPrediction}.
#'
#' @family Flux analysis tools
#' @export
setMethod("show", signature(object = "FluxPrediction"),
          function(object) {
            cat("algorithm:             ", object@algorithm, "\n")
            cat("generic status:        ", object@stat_term, "\n")
            cat("solver status message: ", object@ok_term, "\n")
            cat("Objective fct. value:  ", object@obj, "\n")
            cat("Secondary objective:   ", object@obj_sec, "\n")
          }
)



