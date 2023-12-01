#' Check reaction IDs and Indices
#'
#' Checks whether reaction IDs or indices are part of (valid) for a specific
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

  # Indices provided
  if(is.numeric(react)) {
    checkRes <- react <= react_num(model)
  }

  # IDs provided
  if(is.character(react)) {
    checkRes <- react %in% model@react_id
  }

  return(checkRes)
}


#' Check metabolite IDs and Indices
#'
#' Checks whether metabolite IDs or indices are part of (valid) for a specific
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

  # Indices provided
  if(is.numeric(met)) {
    checkRes <- met <= met_num(model)
  }

  # IDs provided
  if(is.character(met)) {
    checkRes <- met %in% model@met_id
  }

  return(checkRes)
}


#' Check gene IDs and Indices
#'
#' Checks whether gene IDs or indices are part of (valid) for a specific
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

  # Indices provided
  if(is.numeric(gene)) {
    checkRes <- gene <= gene_num(model)
  }

  # IDs provided
  if(is.character(gene)) {
    checkRes <- gene %in% model@allGenes
  }

  return(checkRes)
}

#' Check compartment IDs and Indices
#'
#' Checks whether compartment IDs or indices are part of (valid) for a specific
#' model.
#'
#' @param model Model of class \link{modelorg}
#' @param comp A character vector specifying the compartment IDs or a integer
#' vector providing the compartment indices in the model.
#'
#' @return A logical vector; TRUE if ID/index is valid, FALSE otherwise.
#'
#' @export
checkCompartmentId <- function(model, comp) {
  if (!is(model, "modelorg")) {
    stop("Argument 'model' needs an model of class 'modelorg'.")
  }

  if(!(is.numeric(comp) && all(comp %% 1 == 0) && all(comp > 0) || is.character(comp))) {
    stop("Argument 'comp' needs to be either a non-negative integer vector or a character vector")
  }

  checkRes <- rep(FALSE, length(comp))

  # Indices provided
  if(is.numeric(comp)) {
    checkRes <- comp <= comp_num(model)
  }

  # IDs provided
  if(is.character(comp)) {
    checkRes <- comp %in% model@mod_compart
  }

  return(checkRes)
}

#' Check subsystem IDs and Indices
#'
#' Checks whether subsystem IDs or indices are part of (valid) for a specific
#' model.
#'
#' @param model Model of class \link{modelorg}
#' @param subsystem A character vector specifying the subsystem IDs or a integer
#' vector providing the subsystem indices in the model.
#'
#' @return A logical vector; TRUE if ID/index is valid, FALSE otherwise.
#'
#' @export
checkSubsystemId <- function(model, subsystem) {
  if (!is(model, "modelorg")) {
    stop("Argument 'model' needs an model of class 'modelorg'.")
  }

  if(!(is.numeric(subsystem) && all(subsystem %% 1 == 0) && all(subsystem > 0) || is.character(subsystem))) {
    stop("Argument 'subsystem' needs to be either a non-negative integer vector or a character vector")
  }

  checkRes <- rep(FALSE, length(subsystem))

  # Indices provided
  if(is.numeric(subsystem)) {
    checkRes <- subsystem <= comp_num(model)
  }

  # IDs provided
  if(is.character(subsystem)) {
    checkRes <- subsystem %in% model@subSys_id
  }

  return(checkRes)
}
