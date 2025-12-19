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
#' The idea of flux coupling in flux balance models of multi-species communities
#' was first introduced by Heinken et al. (2013). The idea is to limit the
#' reactions fluxes in order to prevent that a reaction in one organism is
#' solely used (i.e., carries a non-zero flux) to benefit another organism in
#' the community and not the organism that produced the enzyme for the specific
#' reaction. Therefore, new linear constraints are added to the flux balance
#' model where the absolute reaction flux bounds of organism *j* is proportional
#' to the biomass formation of organism *j*. The coupling constraints are
#' defined as:
#'
#' \deqn{-c v_{b,j} - u \leq v_{i,j} \leq c v_{b,j} + u}
#'
#' where \eqn{v_{i,j}} is the flux through reaction *i* in organism *j*,
#' \eqn{v_{b,j}} the biomass reaction of organism *j*. *c* and *u* are
#' parameters for the coupling constraints that define intercept (*u*) and slope
#' (*c*) (see figure).
#'
#' \if{html}{\figure{coupling_constraints.svg}{options: width=300 alt="Coupling Constraints"}}
#' \if{latex}{\figure{coupling_constraints.pdf}{options: width=1.5in}}
#'
#' @references
#' - Heinken A, Sahoo S, Fleming RMT, Thiele I. Systems-level characterization
#' of a host-microbe metabolic symbiosis in the mammalian gut. Vol. 4, Gut
#' Microbes; 2013.
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
