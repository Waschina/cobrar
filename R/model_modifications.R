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
#' @returns An updated model of class \link{modelorg}
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

#' Remove reactions from a model
#'
#' This function removes specified reactions from a model.
#'
#' @param model Model of class \link{modelorg}
#' @param react A character vector stating the reaction IDs in a model or a
#' numeric vector providing the reaction indexes.
#' @param rm_met Logical. Should metabolites, which are singletons after the
#' reaction removal, be deleted as well?
#'
#' @returns An updated model of class \link{modelorg}
#'
#' @export
rmReact <- function(model, react, rm_met = TRUE) {
  if(length(react) == 0)
    return(model)

  if(!all(checkReactId(model, react))) {
    stop("Please check your reaction IDs/indexes in argument 'react'.")
  }

  react <- react_pos(model, react)

  model@S        <- model@S[,-react, drop = FALSE]
  model@obj_coef <- model@obj_coef[-react]
  model@subSys   <- model@subSys[-react,, drop = FALSE]

  model@react_id   <- model@react_id[-react]
  model@react_name <- model@react_name[-react]
  model@react_comp <- model@react_comp[-react]
  model@lowbnd     <- model@lowbnd[-react]
  model@uppbnd     <- model@uppbnd[-react]
  model@react_attr <- model@react_attr[-react,, drop = FALSE]

  model@gprRules <- model@gprRules[-react]
  model@genes    <- model@genes[-react]

  if(rm_met) {
    TMPmat <- as(model@S[,-react], "TsparseMatrix")
    metrm <- which(!(1:met_num(model) %in% (TMPmat@i+1)))
    if(length(metrm) > 0) {

      model@S        <- model@S[-metrm,, drop = FALSE]
      model@met_id   <- model@met_id[-metrm]
      model@met_name <- model@met_name[-metrm]
      model@met_comp <- model@met_comp[-metrm]
      model@met_attr <- model@met_attr[-metrm,, drop = FALSE]

    }
  }

  return(model)
}

#' Remove genes from a model
#'
#' This function removes specified genes from a model, and optionally also
#' reactions and metabolites inaccessible after gene knock outs.
#'
#' @param model Model of class \link{modelorg}
#' @param gene A character vector stating the reaction IDs in a model or a
#' numeric vector providing the reaction indexes.
#' @param rm_react Logical. Should reaction, which are inaccessible after the
#' gene knock outs, be deleted as well?
#' @param rm_met Logical. Should metabolites, which are singletons after the
#' reaction removal, be deleted as well?
#'
#' @returns An updated model of class \link{modelorg}
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#' mod
#'
#' # create a double gene knock-out mutant
#' mod_KO <- rmGene(mod, c("b4152","b0116"))
#' mod_KO
#'
#'
#' @export
rmGene <- function(model, gene, rm_react = TRUE, rm_met = TRUE) {
  if(length(gene) == 0)
    return(model)

  if(!all(checkGeneId(model, gene))) {
    stop("Please check your gene IDs/indexes in argument 'gene'.")
  }

  gene <- gene_pos(model, gene)

  rmReactions <- geneDel(model, gene)

  # rm gene parts
  model@allGenes <- model@allGenes[-gene]
  model@allGenes_name <- model@allGenes_name[-gene]
  model@genes <- lapply(model@genes,
                        function(x) ifelse(x %in% gene, NA_character_,x))

  # rm reaction (and metabolites)
  if(rm_react)
    model <- rmReact(model, rmReactions)

  return(model)
}
