#------------------------------------------------------------------------------#
# Functions to modify the structure / properties of a modelorg                 #
#------------------------------------------------------------------------------#

#' Change flux bounds
#'
#' The function changes either upper bounds, lower bounds, or both for specific
#' reactions.
#'
#' @param model Model of class \link{modelorg}
#' @param react A character vector stating the reaction IDs in a model or a
#' numeric vector providing the reaction indexes.
#' @param lb A numeric vector giving the new lower flux bounds for reactions
#' \code{react}. If \code{lb} is of length 1, the same value will be used for
#' all reactions.
#' @param ub A numeric vector giving the new upper flux bounds for reactions
#' \code{react}. If \code{ub} is of length 1, the same value will be used for
#' all reactions.
#'
#' @export
changeBounds <- function(model, react, lb = NULL, ub = NULL) {

  stopifnot("lb must be numeric" = is.null(lb) || is.numeric(lb),
            "ub must be numeric" = is.null(ub) || is.numeric(ub))

  if(!all(checkReactId(model, react))) {
    stop("Please check your reaction IDs/indexes in argument 'react'.")
  }

  # the actual change
  if(!is.null(lb)) {
    model@lowbnd[match(react, model@react_id)] <- lb
  }
  if(!is.null(ub)) {
    model@lowbnd[match(react, model@react_id)] <- ub
  }

  return(model)
}
