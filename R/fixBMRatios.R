#' Fixate biomass ratios
#'
#' Fixate the biomass ratios in a community models according to models' relative
#' abundances.
#'
#' @param model Community model of class \link{ModelComm}
#' @param BMreact IDs of individual biomass reactions
#' @param tolerance Tolerated deviation of ratios
#'
#' @description
#' Biomass ratios are fixed via introducing new linear constraints.
#'
#' @export
fixBMRatios <- function(model, BMreact = guessBMReaction(model), tolerance = 0) {
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
  BMreact <- BMreact[match(bmIdx,model@community$index)]

  n <- length(BMreact)

  model <- addConstraint(model,
                         react = lapply(1:n, function(x) BMreact),
                         coeff = lapply(1:n, function(x) {
                           cidx <- rep(model@community$abun[x], n)
                           cidx[x] <- cidx[x]-1
                           return(cidx)
                         }),
                         rtype = rep(ifelse(tolerance==0,"E","D"),n),
                         lb = rep(-tolerance, n), ub = rep(tolerance, n))

  return(model)
}
