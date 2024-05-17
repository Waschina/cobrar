#' Parsimonious Flux Balance Analysis (pFBA)
#'
#' Performs parsimonious FBA as describe by Lewis et al. 2010.
#'
#' @param model Model of class \link{ModelOrg}
#' @param costcoeffw,costcoefbw A numeric vector containing cost coefficients
#' for all variables/reactions (forward direction: 'costcoeffw'; backward
#' direction: 'costcoefbw'). If set to NULL, all cost coefficients are set to 1,
#' so that all variables have the same impact on the objective function.
#'
#' @returns
#' Returned reduced costs ('redCosts') correspond to the optimization of the
#' initial linear program (LP), which is basically the initial FBA to calculate
#' the optimal value of the objective function that is defined in 'model'.
#'
#' @references
#' N. E. Lewis et al., “Omic data from evolved E. coli are consistent with
#' computed optimal growth from genome‐scale models,” Molecular Systems Biology,
#' vol. 6, no. 1. EMBO, Jan. 2010. doi: 10.1038/msb.2010.47.
#'
#' @seealso [pfbaHeuristic()]
#'
#' @family Flux prediction algorithms
#' @export
pfba <- function(model, costcoeffw = NULL, costcoefbw = NULL) {

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

  # optimization
  lp_ok   <- solveLp(LPprob)
  lp_stat <- getSolStat(LPprob)

  # get solution (objective value)
  objRes <- getObjValue(LPprob)
  redCosts <- getRedCosts(LPprob)

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

  # add new row to maintain optimal objective function value
  newRow <- Matrix(c(model@obj_coef, -model@obj_coef), nrow = 1, sparse = TRUE)
  newS <- rbind(newS,newRow)

  # new user constraint matrix
  newConstrMat <- cbind(model@constraints@coeff, -model@constraints@coeff)

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

  if(length(bw) != length(fw) || length(fw) != nc)
    stop("The weight vectors 'costcoeffw' and 'costcoefbw' must both be of length equal to the number of reactions.")

  LPprobNew <- new(paste0("LPproblem_",COBRAR_SETTINGS("SOLVER")),
                   name = paste0("LP_", model@mod_id,"_pFBA"),
                   method = COBRAR_SETTINGS("METHOD"))

  loadLPprob(LPprobNew,
             nCols = ncol(newS),
             nRows = nr+1+constraint_num(model),
             mat   = rbind(newS,
                           newConstrMat),
             ub    = newUB,
             lb    = newLB,
             obj   = c(fw,bw),
             rlb   = c(rep(0, nr), objRes,
                       model@constraints@lb),
             rtype = c(rep("E", nr+1),
                       model@constraints@rtype),
             lpdir = "min",
             rub   = c(rep(NA, met_num(model)+1),
                       model@constraints@ub),
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

  return(new("FluxPrediction",
             algorithm = "pFBA",
             ok = lp_ok$code,
             ok_term = lp_ok$term,
             stat = lp_stat$code,
             stat_term = lp_stat$term,
             obj = objRes,
             obj_sec = objResMTF,
             fluxes = lp_fluxes,
             redCosts = redCosts))

}
