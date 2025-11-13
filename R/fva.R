#' Flux Variability Analysis (FVA)
#'
#' Perform Flux Variability Analysis with or without relaxed optimality
#' constraint
#'
#' @param model Model of class \link{ModelOrg}
#' @param react Character vector of reaction IDs tested for flux variability. If
#' NULL, all reactions are tested.
#' @param opt.factor Numeric value > 0 to define the required fraction of the
#' objective function value. E.g. 0.8 sets the constraint, that in the flux
#' variability analysis, the objective function value must at least be 80% of
#' the original optimal value.
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#'
#' # Get flux variability for all exchange reactions
#' fvares <- fva(mod, react = mod@react_id[grepl("^EX_",mod@react_id)],
#'               opt.factor = 0.9)
#' fvares
#'
#' @family Flux prediction algorithms
#' @export
fva <- function(model, react = NULL, opt.factor = 1) {

  if(is.null(react) || length(react) == 0)
    react <- 1:react_num(model)

  if(!all(checkReactId(model, react))) {
    stop("Please check your reaction IDs/indices in argument 'react'.")
  }
  
  if(all(model@obj_coef == 0))
    warning("No objective function defined in the model.")

  #----------------------------------------------------------------------------#
  # First: Basic FBA to find objective value                                   #
  #----------------------------------------------------------------------------#

  # Init LP
  LPprob <- new(paste0("LPproblem_",COBRAR_SETTINGS("SOLVER")),
                name = paste0("LP_", model@mod_id),
                method = COBRAR_SETTINGS("METHOD"))

  loadLPprob(LPprob,
             nCols = react_num(model),
             nRows = met_num(model)+constraint_num(model),
             mat   = rbind(model@S, model@constraints@coeff),
             ub    = ifelse(abs(model@uppbnd)>COBRAR_SETTINGS("MAXIMUM"),
                            sign(model@uppbnd)*COBRAR_SETTINGS("MAXIMUM"),
                            model@uppbnd),
             lb    = ifelse(abs(model@lowbnd)>COBRAR_SETTINGS("MAXIMUM"),
                            sign(model@lowbnd)*COBRAR_SETTINGS("MAXIMUM"),
                            model@lowbnd),
             obj   = model@obj_coef,
             rlb   = c(rep(0, met_num(model)),
                       model@constraints@lb),
             rtype = c(rep("E", met_num(model)),
                       model@constraints@rtype),
             lpdir = substr(model@obj_dir,1,3),
             rub   = c(rep(NA, met_num(model)),
                       model@constraints@ub),
             ctype = NULL
  )

  # optimization
  lp_ok   <- solveLp(LPprob)
  lp_stat <- getSolStat(LPprob)

  # get solution (objective value)
  objRes <- getObjValue(LPprob)

  if(is.na(objRes)) {
    stop(paste("There is no feasible FBA solution and no flux variablity can be calculated. Solver status:", lp_stat$term))
  }

  #----------------------------------------------------------------------------#
  # Preparation for FVA: fix range of value of objective function as new       #
  # constraint (=row)                                                          #
  #----------------------------------------------------------------------------#
  addSingleConstraint(LPprob,
                      model@obj_coef,
                      min(objRes*opt.factor[1], objRes),
                      max(objRes*opt.factor[1], objRes),
                      ifelse(opt.factor[1] == 1, "E","D"))


  react <- react_pos(model, react)

  res <- list()
  for(i in 1:length(opt.factor)) {
    if(i > 1) {
      setRowsBnds(LPprob,
                  i = met_num(model)+constraint_num(model)+1,
                  objRes*opt.factor[i],
                  objRes,
                  ifelse(opt.factor[i] == 1, "E","D"))
    }

    res[[i]] <- fvaJob(LPprob, react)
    res[[i]] <- cbind(data.frame(react = model@react_id[react],
                                 growth.fraction = opt.factor[i]),
                      res[[i]])
  }

  res <- do.call("rbind",res)

  #----------------------------------------------------------------------------#
  # Delete LP-Problem and free associated memory                               #
  #----------------------------------------------------------------------------#
  deleteLP(LPprob)

  return(res)
}
