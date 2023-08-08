#' @export
fva <- function(model, react = NULL, opt.factor = 1) {

  if(is.null(react) || length(react) == 0)
    react <- 1:react_num(model)

  if(!all(checkReactId(model, react))) {
    stop("Please check your reaction IDs/indexes in argument 'react'.")
  }

  #----------------------------------------------------------------------------#
  # First: Basic FBA to find objective value                                   #
  #----------------------------------------------------------------------------#

  # Init LP
  LPprob <- new(paste0("LPproblem_",COBRAR_SETTINGS("SOLVER")),
                name = paste0("LP_", model@mod_id),
                method = COBRAR_SETTINGS("METHOD"))

  loadLPprob(LPprob,
             nCols = react_num(model),
             nRows = met_num(model),
             mat   = model@S,
             ub    = model@uppbnd,
             lb    = model@lowbnd,
             obj   = model@obj_coef,
             rlb   = rep(0, met_num(model)),
             rtype = rep("E", met_num(model)),
             lpdir = COBRAR_SETTINGS("OPT_DIRECTION"),
             rub   = NULL,
             ctype = NULL
  )

  # optimization
  lp_ok   <- solveLp(LPprob)
  lp_stat <- getSolStat(LPprob)

  # get solution (objective value)
  objRes <- getObjValue(LPprob)

  #----------------------------------------------------------------------------#
  # Preparation for FVA: fix range of value of objective function as new       #
  # constraint (=row)                                                          #
  #----------------------------------------------------------------------------#
  addSingleConstraint(LPprob,
                      model@obj_coef,
                      objRes*opt.factor[1],
                      objRes,
                      ifelse(opt.factor[1] == 1, "E","D"))

  react <- react_pos(model, react)

  res <- list()
  for(i in 1:length(opt.factor)) {
    if(i > 1) {
      setRowsBnds(LPprob,
                  i = met_num(model)+1,
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
  return(res)

}
