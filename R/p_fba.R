#' @export
p_fba <- function(model, costcoeffw = NULL, costcoefbw = NULL) {

  if(!is.null(costcoeffw) && !is.numeric(costcoeffw))
    stop("Argument 'costcoeffw' must be a numeric vector")
  if(!is.null(costcoefbw) && !is.numeric(costcoefbw))
    stop("Argument 'costcoefbw' must be a numeric vector")

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
  # New LP for minimalization of total flux                                    #
  #----------------------------------------------------------------------------#

  nc <- react_num(model)
  nr <- met_num(model)

  # Now make all reactions irreversible and operate solely from LHS->RHS
  aS  <- model@S
  aLB <- ifelse(model@lowbnd > 0, model@lowbnd, 0)
  aUB <- ifelse(model@uppbnd > 0, model@uppbnd, 0)

  bS <- -model@S
  bLB <- ifelse(-model@uppbnd > 0, -model@uppbnd, 0)
  bUB <- ifelse(-model@lowbnd > 0, -model@lowbnd, 0)

  newS  <- cbind(aS, bS)
  newLB <- c(aLB,bLB)
  newUB <- c(aUB,bUB)

  # add new roe to maintain optimal objective function value
  newRow <- Matrix(c(model@obj_coef, -model@obj_coef), nrow = 1, sparse = TRUE)
  newS <- rbind(newS,newRow)

  # Init new LP
  if(is.null(costcoeffw)) {
    fw <- rep(1, nc)
  } else {
    fw <- costcoeffw
  }

  if(is.null(costcoefbw)) {
    bw <- fw
  } else {
    bw <- costcoefbw
  }

  LPprobNew <- new(paste0("LPproblem_",COBRAR_SETTINGS("SOLVER")),
                   name = paste0("LP_", model@mod_id,"_pFBA"),
                   method = COBRAR_SETTINGS("METHOD"))

  loadLPprob(LPprobNew,
             nCols = ncol(newS),
             nRows = nr+1,
             mat   = newS,
             ub    = newUB,
             lb    = newLB,
             obj   = c(fw,bw),
             rlb   = c(rep(0, nr),objRes),
             rtype = rep("E", nr+1),
             lpdir = "min",
             rub   = NULL,
             ctype = NULL
  )

  # optimization
  lp_ok   <- solveLp(LPprobNew)
  lp_stat <- getSolStat(LPprobNew)

  # get solution (objective value)
  objResMTF <- getObjValue(LPprobNew)
  lp_fluxes <- getColsPrimal(LPprobNew)
  fwflx <- lp_fluxes[1:nc]; bwflx <- lp_fluxes[(nc+1):(nc*2)]
  lp_fluxes <- ifelse(bwflx > fwflx, -bwflx, fwflx)

  return(list(ok = lp_ok$code,
              ok_term = lp_ok$term,
              stat = lp_stat$code,
              stat_term = lp_stat$term,
              obj_fba = objRes,
              obj_mtf = objResMTF,
              fluxes = lp_fluxes
  ))

}
