#------------------------------------------------------------------------------#
# Functions to do basic analysis of predicted flux distributions               #
#------------------------------------------------------------------------------#

#' Get metabolite exchange rates
#'
#' Summarize the predicted fluxes for exchange reactions.
#'
#' @param model Model of class \link{ModelOrg}
#' @param sol Flux prediction results as an object of class
#' \link{FluxPrediction}.
#'
#' @returns A data.frame that summarizes the predicted metabolite exchange rates
#' (=fluxes of exchange reactions).
#'
#' @export
getExchanges <- function(model, sol) {

  exind <- findExchReact(model)$react_pos

  df <- data.frame(ID = model@react_id[exind],
                   name = model@react_name[exind],
                   flux = sol@fluxes[exind])

  return(df)
}
