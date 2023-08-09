#' Heuristic parsimonious Flux Balance Analysis (pFBA)
#'
#' Performs a heuristic version of the parsimonious FBA algorithm. See
#' details.
#'
#' @param model Model of class \link{modelorg}
#' @param costcoeffw,costcoefbw A numeric vector containing cost coefficients
#' for all variables/reactions (forward direction: 'costcoeffw'; backward
#' direction: 'costcoefbw'). If set to NULL, all cost coefficients are set to 1,
#' so that all variables have the same impact on the objective function.
#' @param pFBAcoeff Numeric value to weight the minimization of total flux
#' within the combined objective function. See details.
#'
#' @details
#' The exact-solution pFBA algorithm described by Lewis et al. 2010 consists of
#' two optimization steps: (1) A basic flux balance analysis is performed to
#' obtain the optimal value of the objective function (e.g., growth rate). (2)
#' The objective function from (1) is transformed into a constraint where the
#' value of the function is fixed to the determined optimal value. A new
#' objective is defined that minimizes the absolute sum of fluxes through the
#' metabolic network.\cr
#' The here implemented heuristic pFBA performs only one optimization step.
#' Therefore, the original objective function is combined with a term that
#' minimizes the absolute (weighted) sum of fluxes:
#' \deqn{max: \sum_{i=1}^{n}(c_i v_i) - p \sum_{i=1}^{n}(w_i |v_i|)}
#' for maximization direction or
#' \deqn{min: \sum_{i=1}^{n}(c_i v_i) + p \sum_{i=1}^{n}(w_i |v_i|)}
#' if the original objective function is minimized.\cr
#' Here,
#' \eqn{c_i} is the original objective coefficient of reaction \eqn{i},
#' \eqn{v_i} the flux through reaction \eqn{i},
#' \eqn{n} the number of reactions in the model, and
#' \eqn{w_i} the weight of reaction \eqn{i} (see arguments 'costcoeffw',
#' 'costcoefbw').
#' The scalar parameter \eqn{p} (argument 'pFBAcoeff') defines the weighting of
#' the total flux minimization. Increasing this value increases the weight of
#' total flux minimization, possibly at costs for the value of the objective
#' function defined in 'model' (e.g., flux through biomass reaction).\cr
#' This heuristic implementation of a pFBA is the core of the gap-filling
#' algorithm of 'gapseq' (Zimmermann et al. 2021).
#'
#'
#' @references
#' N. E. Lewis et al., “Omic data from evolved E. coli are consistent with
#' computed optimal growth from genome‐scale models,” Molecular Systems Biology,
#' vol. 6, no. 1. EMBO, Jan. 2010. doi: 10.1038/msb.2010.47.
#'
#' Zimmermann, C. Kaleta, and S. Waschina, “gapseq: informed prediction of
#' bacterial metabolic pathways and reconstruction of accurate metabolic
#' models,” Genome Biology, vol. 22, no. 1. Springer Science and Business Media
#' LLC, Mar. 10, 2021. doi: 10.1186/s13059-021-02295-1.
#'
#' @export
pfba_heuristic <- function(model, costcoeffw = NULL, costcoefbw = NULL,
                            pFBAcoeff = 1e-6) {

  if(!is.null(costcoeffw) && !is.numeric(costcoeffw))
    stop("Argument 'costcoeffw' must be a numeric vector")
  if(!is.null(costcoefbw) && !is.numeric(costcoefbw))
    stop("Argument 'costcoefbw' must be a numeric vector")


  nc <- react_num(model)
  nr <- met_num(model)

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

  # define new objective function (combined original objevtive with factorized mtf)
  if(COBRAR_SETTINGS("OPT_DIRECTION") == "max")
    newObj <- c(model@obj_coef, -model@obj_coef) - pFBAcoeff * c(fw, bw)
  if(COBRAR_SETTINGS("OPT_DIRECTION") == "min")
    newObj <- c(model@obj_coef, -model@obj_coef) + pFBAcoeff * c(fw, bw)

  # define and parametrize LP problem
  LPprob <- new(paste0("LPproblem_",COBRAR_SETTINGS("SOLVER")),
                name = paste0("LP_", model@mod_id,"_pFBA_heuristic"),
                method = COBRAR_SETTINGS("METHOD"))

  loadLPprob(LPprob,
             nCols = ncol(newS),
             nRows = nr,
             mat   = newS,
             ub    = newUB,
             lb    = newLB,
             obj   = newObj,
             rlb   = rep(0, nr),
             rtype = rep("E", nr),
             lpdir = COBRAR_SETTINGS("OPT_DIRECTION")
  )

  # optimization
  lp_ok   <- solveLp(LPprob)
  lp_stat <- getSolStat(LPprob)

  # get solution (objective value)
  objResMTF <- getObjValue(LPprob)

  lp_fluxes <- getColsPrimal(LPprob)
  fwflx <- lp_fluxes[1:nc]; bwflx <- lp_fluxes[(nc+1):(nc*2)]
  lp_fluxes <- ifelse(bwflx > fwflx, -bwflx, fwflx)

  objRes    <- as.numeric(lp_fluxes %*% model@obj_coef)
  redCosts  <- getRedCosts(LPprob)


  return(new("FluxPrediction",
             algorithm = "pFBA (heuristic)",
             ok = lp_ok$code,
             ok_term = lp_ok$term,
             stat = lp_stat$code,
             stat_term = lp_stat$term,
             obj = objRes,
             obj_sec = objResMTF,
             fluxes = lp_fluxes,
             redCosts = redCosts[1:nc]))

}
