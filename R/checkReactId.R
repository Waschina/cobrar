#' Check reaction IDs and Indexes
#'
#' Checks whether reaction IDs or indexes are part of (valid) for a specific
#' model.
#'
#' @param object Model of class \link{modelorg}
#' @param react A character vector specifying the reaction IDs or a integer
#' vector providing the reaction indices in the model.
#'
#' @return A logical vector; TRUE if ID/index is valid, FALSE otherwise.
#'
#' @export
checkReactId <- function(object, react) {
  if (!is(object, "modelorg")) {
    stop("Argument 'object' needs an object of class 'modelorg'.")
  }

  if(!(is.numeric(react) && all(react %% 1 == 0) && all(react > 0) || is.character(react))) {
    stop("Argument 'react' needs to be either a non-negative integer vector or a character vector")
  }

  checkRes <- rep(FALSE, length(react))

  # Indexes provided
  if(is.numeric(react)) {
    checkRes <- react <= react_num(object)
  }

  # IDs provided
  if(is.character(react)) {
    checkRes <- react %in% object@react_id
  }

  return(checkRes)
}
