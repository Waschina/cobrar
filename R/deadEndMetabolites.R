#' Identify dead end metabolites
#'
#' Searches a metabolic network for metabolites that can be produced but not
#' consumed, and vice versa.
#'
#' @param object Model of class \link{modelorg}
#'
#' @returns A list with two elements: "dem" is a character vector with the IDs
#' of identified dead end metabolites. "der" is a character vector with IDs of
#' reactions, that involve at least one of the dead end metabolites.
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#' deadEndMetabolites(mod)
#'
#' @note
#' The algorithm is adapted from the original 'sybil' package, which in turn is
#' adapted from the implementation in 'cobratoolbox'. Important detail: The
#' algorithm function in the 'sybil' package may return different results than
#' this cobrar implementation. This is because the 'sybil' function does not
#' take into account the possibility that an irreversible reaction is defined in
#' the direction RHS to LHS (lower bound > 0 and upper bound = 0).
#'
#' @export
deadEndMetabolites <- function(object) {

  St <- object@S
  LBt <- object@lowbnd
  UBt <- object@uppbnd
  St[which(abs(St) < COBRAR_SETTINGS("TOLERANCE"), arr.ind = TRUE)] <- 0

  # turn reactions, which are only operating from RHS to LHS
  # this simplifies the algorithm below
  ind_turn <- which(LBt < 0 & UBt <= 0)
  if(length(ind_turn) > 0) {
    St[,ind_turn] <- -St[,ind_turn]
    tmp <- LBt
    LBt[ind_turn] <- -UBt[ind_turn]
    UBt[ind_turn] <- -tmp[ind_turn]
    rm(tmp)
  }

  dem <- rep(FALSE, met_num(object))
  der <- rep(FALSE, react_num(object))
  keep_going <- TRUE
  while(keep_going) {
    newDem <- apply(St, 1, function(x) {
      NZcol <- x != 0
      !(sum(NZcol)>=2 && (
        (any(x>0) && any(x<0)) ||
          any(LBt[NZcol] < 0))
        )
    })
    newDem[dem] <- FALSE
    dem[newDem] <- TRUE

    newDer <- apply(St[newDem,, drop = FALSE],2,function(x) any(x != 0))
    der[newDer] <- TRUE
    if(any(newDer)) {
      St[,newDer, drop = FALSE] <- 0
    }
    # print(sum(dem))
    if(sum(newDem) == 0)
      keep_going <- FALSE
  }
  dem <- object@met_id[dem]
  der <- object@react_id[der]

  return(list(dem = dem, der = der))
}
