#' Couple reaction flux bounds to biomass production
#'
#' Resets the absolute flux bounds of reactions proportional to the flux through
#' the biomass reaction.
#'
#' @param model Community model of class \link{ModelComm}
#' @param BMreact IDs of individual biomass reactions
#' @param cpl_c Coupling constraint \code{c}
#' @param cpl_u Coupling constraint \code{u}
#'
#' @details
#' See (and cite) Heinken et al. (2013) Gut Microbes (doi: 10.4161/gmic.22370).
#'
#' @export
fluxBMCoupling <- function(model, BMreact = guessBMReaction(model),
                           cpl_c = 400, cpl_u = 0.01) {
  if(length(BMreact) != nrow(model@community))
    stop("Length of 'BMreact' not equal to number of organisms in the community.")

  bmIdx <- unlist(regmatches(BMreact,gregexpr("^M[0-9]+_", BMreact)))
  bmIdx <- gsub("_$","",bmIdx)

  if(any(duplicated(bmIdx)))
    stop("There are two or more biomass reactions for at least on organism, which is not allowed. Please check 'BMreact'.")

  if(!all(bmIdx %in% model@community$index))
    stop("At least one provided biomass reaction ID cannot be matched to a community organism.")

  # make sure biomass reactions are in the same order as organisms within the
  # community
  names(BMreact) <- bmIdx
  BMreact <- BMreact[match(bmIdx,model@community$index)]
  n <- length(BMreact)

  # identify reactions to be coupled (exclude organisms' exchange reactions)
  mrxns <- model@react_id[grepl("^M[0-9]+_", model@react_id)]
  mrxns <- mrxns[!grepl("^M[0-9]+_EX_", mrxns)]
  mrxns <- mrxns[!(mrxns %in% BMreact)]
  mrxnsIdx <- unlist(regmatches(mrxns,gregexpr("^M[0-9]+_", mrxns)))
  mrxnsIdx <- gsub("_$","",mrxnsIdx)
  reaBMpair <- mapply(function(x,y) {c(x,y)}, BMreact[mrxnsIdx], mrxns,
                      SIMPLIFY = FALSE)

  nr <- length(mrxns)
  reaCoefL <- mapply(function(x,y) {c(x,y)}, rep(cpl_c,nr), rep(-1,nr),
                     SIMPLIFY = FALSE)
  reaCoefU <- mapply(function(x,y) {c(x,y)}, rep(-cpl_c,nr), rep(-1,nr),
                     SIMPLIFY = FALSE)

  model <- addConstraint(model,
                         react = c(reaBMpair, reaBMpair),
                         coeff = c(reaCoefL, reaCoefU),
                         rtype = c(rep("L",nr),rep("U",nr)),
                         lb = rep(-cpl_u,nr*2),
                         ub = rep(cpl_u,nr*2))
  return(model)
}
