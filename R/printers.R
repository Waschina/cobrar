#' Print reaction(s)
#'
#' Print the equations of reactions.
#'
#' @param model Model of class \link{ModelOrg}
#' @param react A character vector specifying the reaction IDs or a integer
#' vector providing the reaction indices in the model.
#' @param use.ids Boolean. Indicating whether metabolite IDs should be printed
#' instead of metabolite names.
#'
#' @return A character vector with the individual reaction equations.
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#' # print reaction specified by index
#' printReaction(mod, react = 2)
#' # print reaction specified by ID
#' printReaction(mod, react = "PFL")
#' # print reaction with metabolite IDs instead of metabolite names
#' printReaction(mod, react = "PFL", use.ids = TRUE)
#' # print multiple reactions at once
#' printReaction(mod, react = c(2,8))
#'
#' @family Model characteristics
#' @export
printReaction <- function(model, react, use.ids = FALSE) {
  check <- checkReactId(model, react = react)
  if(any(!check)) {
    stop("check argument react")
  }

  cind <- react_pos(model, react)
  mat <- model@S[, cind, drop = FALSE]
  nnz <- apply(mat, 2, "!=", 0)
  reaction <- character(length(cind))

  for (j in seq(along = cind)) {

    if(use.ids) {
      met <- model@met_id[nnz[, j]]
    } else {
      met <- model@met_name[nnz[, j]]
    }

    nzv <- mat[, j][nnz[, j]]

    ed <- nzv < 0
    pd <- nzv > 0

    if (sum(ed) > 0) {
      educt   <- paste(paste("(", abs(nzv[ed]), ")", sep = ""),
                       met[ed], collapse = " + ")
    }
    else {
      educt = ""
    }

    if (sum(pd) > 0) {
      product <- paste(paste("(", nzv[pd], ")", sep = ""),
                       met[pd], collapse = " + ")
    }
    else {
      product = ""
    }

    arrow <- " <==> "
    if(model@lowbnd[cind[j]] >= 0 & model@uppbnd[cind[j]] > 0)
      arrow <- " --> "
    if(model@lowbnd[cind[j]] < 0 & model@uppbnd[cind[j]] <= 0)
      arrow <- " <-- "

    reaction[j] <- paste(educt, product, sep = arrow)
  }

  # names(reaction) <- react

  return(reaction)
}

#' Print Constraint(s)
#'
#' Generate strings to summarize metabolic model constraints.
#'
#' @param model Model of class \link{ModelOrg}
#' @param ind Integer vector with the indices of the constraints to be printed
#'
#' @family Model characteristics
#' @export
printConstraint <- function(model, ind = NULL) {
  if(is.null(ind) && constraint_num(model) == 0)
    return(character(0))

  if(is.null(ind))
    ind <- 1:constraint_num(model)

  if(constraint_num(model) == 0 || any(!(ind %in% 1:constraint_num(model)))) {
    stop("Invalid index or indices for constraints.")
  }
  res <- sapply(ind, FUN = function(i) constraint2string(model, i))

  return(res)
}
