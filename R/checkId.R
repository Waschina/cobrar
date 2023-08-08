#' Check reaction IDs and Indexes
#'
#' Checks whether reaction IDs or indexes are part of (valid) for a specific
#' model.
#'
#' @param model Model of class \link{modelorg}
#' @param react A character vector specifying the reaction IDs or a integer
#' vector providing the reaction indices in the model.
#'
#' @return A logical vector; TRUE if ID/index is valid, FALSE otherwise.
#'
#' @export
checkReactId <- function(model, react) {
  if (!is(model, "modelorg")) {
    stop("Argument 'model' needs an model of class 'modelorg'.")
  }

  if(!(is.numeric(react) && all(react %% 1 == 0) && all(react > 0) || is.character(react))) {
    stop("Argument 'react' needs to be either a non-negative integer vector or a character vector")
  }

  checkRes <- rep(FALSE, length(react))

  # Indexes provided
  if(is.numeric(react)) {
    checkRes <- react <= react_num(model)
  }

  # IDs provided
  if(is.character(react)) {
    checkRes <- react %in% model@react_id
  }

  return(checkRes)
}


#' Check metabolite IDs and Indexes
#'
#' Checks whether metabolite IDs or indexes are part of (valid) for a specific
#' model.
#'
#' @param model Model of class \link{modelorg}
#' @param react A character vector specifying the metabolite IDs or a integer
#' vector providing the metabolite indices in the model.
#'
#' @return A logical vector; TRUE if ID/index is valid, FALSE otherwise.
#'
#' @export
checkMetId <- function(model, met) {
  if (!is(model, "modelorg")) {
    stop("Argument 'model' needs an model of class 'modelorg'.")
  }

  if(!(is.numeric(met) && all(met %% 1 == 0) && all(met > 0) || is.character(met))) {
    stop("Argument 'met' needs to be either a non-negative integer vector or a character vector")
  }

  checkRes <- rep(FALSE, length(met))

  # Indexes provided
  if(is.numeric(met)) {
    checkRes <- met <= met_num(model)
  }

  # IDs provided
  if(is.character(met)) {
    checkRes <- met %in% model@met_id
  }

  return(checkRes)
}


#' Check gene IDs and Indexes
#'
#' Checks whether gene IDs or indexes are part of (valid) for a specific
#' model.
#'
#' @param model Model of class \link{modelorg}
#' @param gene A character vector specifying the gene IDs or a integer
#' vector providing the gene indices in the model.
#'
#' @return A logical vector; TRUE if ID/index is valid, FALSE otherwise.
#'
#' @export
checkGeneId <- function(model, gene) {
  if (!is(model, "modelorg")) {
    stop("Argument 'model' needs an model of class 'modelorg'.")
  }

  if(!(is.numeric(gene) && all(gene %% 1 == 0) && all(gene > 0) || is.character(gene))) {
    stop("Argument 'gene' needs to be either a non-negative integer vector or a character vector")
  }

  checkRes <- rep(FALSE, length(gene))

  # Indexes provided
  if(is.numeric(gene)) {
    checkRes <- gene <= gene_num(model)
  }

  # IDs provided
  if(is.character(gene)) {
    checkRes <- gene %in% model@allGenes
  }

  return(checkRes)
}
