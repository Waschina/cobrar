#' Flux Balance Analysis
#'
#' Performs basic flux balance analysis (fba)
#'
#' @param model Model of class \link{modelorg}
#'
#' @returns A list with flux predictions (reaction fluxes 'fluxes', reduced costs 'redCosts'), and optimization status ()
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#'
#' # aerobic growth
#' res_aero <- fba(mod)
#' cat(" Growth rate:       ", res_aero$obj,"\n",
#'     "Acetate prodcution:", res_aero$fluxes[mod@react_id == "EX_ac_e"],"\n")
#'
#' mod <- changeBounds(mod, react = "EX_o2_e", lb = 0) # before: -1000
#' res_anaero <- fba(mod)
#' cat(" Growth rate:       ", res_anaero$obj,"\n",
#'     "Acetate prodcution:", res_anaero$fluxes[mod@react_id == "EX_ac_e"],"\n")
#'
#' @export
fba <- function(model) {

  #----------------------------------------------------------------------------#
  # Initializing and defining LP problem                                       #
  #----------------------------------------------------------------------------#
  LPprob <- new(paste0("LPproblem_",COBRAR_SETTINGS("SOLVER")),
                name = paste0("LP_", model@mod_id),
                method = COBRAR_SETTINGS("METHOD"))

  loadLPprob(LPprob,
             nCols = react_num(model),
             nRows = met_num(model)+constraint_num(model),
             mat   = rbind(model@S, model@constraints@coeff),
             ub    = model@uppbnd,
             lb    = model@lowbnd,
             obj   = model@obj_coef,
             rlb   = c(rep(0, met_num(model)),
                       model@constraints@lb),
             rtype = c(rep("E", met_num(model)),
                       model@constraints@rtype),
             lpdir = COBRAR_SETTINGS("OPT_DIRECTION"),
             rub   = c(rep(NA, met_num(model)),
                       model@constraints@ub),
             ctype = NULL
  )

  #----------------------------------------------------------------------------#
  # Optimizing problem                                                         #
  #----------------------------------------------------------------------------#
  lp_ok   <- solveLp(LPprob)
  lp_stat <- getSolStat(LPprob)

  #----------------------------------------------------------------------------#
  # Retrieve predictions                                                       #
  #----------------------------------------------------------------------------#
  objRes <- getObjValue(LPprob)
  lp_fluxes <- getColsPrimal(LPprob)

  redCosts <- getRedCosts(LPprob)


  return(new("FluxPrediction",
             algorithm = "FBA",
             ok = lp_ok$code,
             ok_term = lp_ok$term,
             stat = lp_stat$code,
             stat_term = lp_stat$term,
             obj = objRes,
             obj_sec = NA_real_,
             fluxes = lp_fluxes,
             redCosts = redCosts))
}

