#' Print reaction(s)
#'
#' Print the equations of reactions.
#'
#' @param object Model of class \link{modelorg}
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
#' printReaction(mod, react = "R_PFL")
#' # print reaction with metabolite IDs instead of metabolite names
#' printReaction(mod, react = "R_PFL", use.ids = TRUE)
#' # print multiple reactions at once
#' printReaction(mod, react = c(2,8))
#'
#' @export
printReaction <- function(object, react, use.ids = FALSE) {
  check <- checkReactId(object, react = react)
  if(any(!check)) {
    stop("check argument react")
  }

  cind <- react_pos(object, react)
  mat <- object@S[, cind, drop = FALSE]
  nnz <- apply(mat, 2, "!=", 0)
  reaction <- character(length(cind))

  for (j in seq(along = cind)) {

    if(use.ids) {
      met <- object@met_id[nnz[, j]]
    } else {
      met <- object@met_name[nnz[, j]]
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
    if(object@lowbnd[cind[j]] >= 0 & object@uppbnd[cind[j]] > 0)
      arrow <- " --> "
    if(object@lowbnd[cind[j]] < 0 & object@uppbnd[cind[j]] <= 0)
      arrow <- " <-- "

    reaction[j] <- paste(educt, product, sep = arrow)
  }

  # names(reaction) <- react

  return(reaction)
}
