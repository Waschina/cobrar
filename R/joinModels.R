#' Join multiple metabolic models to form a community
#'
#' Merges multiple metabolic models into one model where each organisms is
#' within an own compartment. Organism compartments are connected by a shared
#' extracellular space (i.e., environment), from where the organisms can take up
#' nutrients an release metabolic by-products. Exchange reactions are attached
#' to the new (shared) extracelluar compartment.
#'
#' @param models A list of objects of class \link{ModelOrg}
#' @param mergeLB,mergeUB Method to use in cases where the individual models
#' differ in the lower/upper bound values for the same exchange reactions. Option
#' "min" uses the minimum LB/UB values found among the individual models for
#' the exchange reactions in the final merged model. "max" acts analogously.
#' "none" will cause the method to stop with an error, if boundaries differ
#' between input models. Further options: "median","mean".
#' @param origExchangeBounds Method to keep/open the flux boundaries of
#' organism-specific exchange reactions. "keep" does not alter the lower or
#' upper bounds of the organism exchange reactions, that connect the respective
#' organism to the shared extracellular space. "open" removes any constraints on
#' organism-specific exchange reactions. "open_lb" and "open_ub" removes only
#' bounds for lower or upper limits, respectively.
#' @param abun Numeric vector to define the relative abundance of the model
#' organisms in the community. Entries of the vector should ideally sum up to 1.
#'
#' @importFrom stats aggregate
#'
#' @export
joinModels <- function(models, mergeLB = "none", mergeUB = "none",
                       origExchangeBounds = "open",
                       abun = rep(1/length(models),length(models))) {

  if(!(mergeLB %in% c("none","min","max","mean","median")))
    stop("Unknown option for 'mergeLB'")

  if(!(mergeUB %in% c("none","min","max","mean","median")))
    stop("Unknown option for 'mergeUB'")

  if(length(abun) != length(models))
    stop("'models' and 'abun' are not of the same length.")

  if(sum(abun) != 1) {
    abun <- abun/sum(abun)
    warning("Values of 'abun' do not sum to 1. Abundances were rescaled.")
  }


  n <- length(models)

  # Initialize community model
  cmod <- new("ModelComm")
  cmod@community <- data.frame(
    index = paste0("M", sprintf(paste0("%0", nchar(n), "d"), 1:n)),
    id = unlist(lapply(models, FUN = function(x) x@mod_id)),
    desc = unlist(lapply(models, FUN = function(x) x@mod_desc)),
    name = unlist(lapply(models, FUN = function(x) x@mod_name)),
    abun = abun
  )

  indeces <- cmod@community$index
  names(models) <- indeces

  #
  # Aggregate lower/upper bounds for new global exchange reactions
  #
  dtex <- lapply(models, function(x) {
    exind <- grep("^EX_", x@react_id)
    metind <- apply(x@S[,exind],2,function(y) which(y != 0))
    res <- data.frame(exrxn = x@react_id[exind],
                      lb = x@lowbnd[exind],
                      ub = x@uppbnd[exind],
                      cpd = x@met_id[metind],
                      cpd.name = x@met_name[metind])
    res <- cbind(res, x@met_attr[metind,])
  })
  dtex <- do.call("rbind", dtex)

  if(mergeLB == "none") {
    nuniq <- aggregate(lb ~ exrxn, data = dtex, function(x) length(unique(x)))
    if(any(nuniq$lb > 1))
      stop(paste("Conflicting lower bounds for exchange reaction:",
                 paste(nuniq[nuniq$lb > 1,"exrxn"], collapse = ", ")))
    mergeLB <- "mean"
  }
  if(mergeUB == "none") {
    nuniq <- aggregate(ub ~ exrxn, data = dtex, function(x) length(unique(x)))
    if(any(nuniq$ub > 1))
      stop(paste("Conflicting upper bounds for exchange reaction:",
                 paste(nuniq[nuniq$ub > 1,"exrxn"], collapse = ", ")))
    mergeUB <- "mean"
  }
  dtexLB <- aggregate(lb ~ exrxn, data = dtex, mergeLB)
  dtexUB <- aggregate(ub ~ exrxn, data = dtex, mergeUB)
  dtexAttr <-
  dtexNew <- merge(dtexLB,dtexUB, by = "exrxn")
  dtexNew <- merge(dtexNew, dtex[!duplicated(dtex$exrxn),
                                 c("exrxn","cpd","cpd.name",
                                   colnames(models[[1]]@met_attr))],
                   by = "exrxn")

  #
  # New stoichiometric matrix
  #
  nrxns <- unlist(lapply(models, react_num))
  rxnoffset <- cumsum(c(0, nrxns[-n]))
  nmets <- unlist(lapply(models, met_num))
  metoffset <- cumsum(c(0, nmets[-n]))

  # `+ nrow(dtexNew)` term for new extracellular metabolites and their global
  # exchange reactions
  cmod@S <- Matrix(0,
                   nrow = sum(nmets) + nrow(dtexNew),
                   ncol = sum(nrxns) + nrow(dtexNew),
                   sparse = TRUE)

  inds <- lapply(models, function(x)  {
    tmp <- which(x@S!=0, arr.ind = TRUE)
    tmp <- cbind(tmp, x@S[tmp])
  })

  for(i in 1:length(inds)) {
    indstmp <- inds[[i]]
    indstmp[,"row"] <- indstmp[,"row"]+metoffset[i]
    indstmp[,"col"] <- indstmp[,"col"]+rxnoffset[i]
    cmod@S[indstmp[,c("row","col")]] <- indstmp[,3]
  }

  #
  # Model compartments
  #
  cmod@mod_compart <- unlist(lapply(indeces,
                                    function(idx) paste0(idx,"_",
                                                         models[[idx]]@mod_compart)))
  cmod@mod_compart_name <- unlist(lapply(indeces,
                                    function(idx) paste0(idx," ",
                                                         models[[idx]]@mod_compart_name)))

  #
  # Metabolites & reaction ids, names, attributes
  #
  cmod@react_id <- unlist(lapply(indeces,
                                 function(idx) paste0(idx,"_",
                                                      models[[idx]]@react_id)))
  cmod@react_name <- unlist(lapply(models, function(x) x@react_name))
  cmod@react_comp <- unlist(lapply(models, function(x) x@react_comp))
  cmod@react_attr <- do.call("rbind",lapply(models, function(x) x@react_attr))

  cmod@met_id <- unlist(lapply(indeces,
                               function(idx) paste0(idx,"_",
                                                    models[[idx]]@met_id)))
  cmod@met_name <- unlist(lapply(models, function(x) x@met_name))
  cmod@met_comp <- unlist(lapply(models, function(x) x@met_comp))
  cmod@met_attr <- do.call("rbind",lapply(models, function(x) x@met_attr))

  #
  # Lower/Upper bounds
  #
  cmod@lowbnd <- unlist(lapply(models, function(x) x@lowbnd))
  cmod@uppbnd <- unlist(lapply(models, function(x) x@uppbnd))

  #
  # Add new global extracellular space
  #
  cmod <- addCompartment(cmod, "e", "shared extracellular space")

  #
  # Add new global exchange reactions
  #
  newExSinds <- matrix(c((1:nrow(dtexNew)) + length(cmod@met_id),
                         (1:nrow(dtexNew)) + length(cmod@react_id)),
                       ncol = 2)
  cmod@S[newExSinds] <- rep(-1, nrow(newExSinds))

  cmod@react_id <- c(cmod@react_id, dtexNew$exrxn)
  cmod@react_name <- c(cmod@react_name,paste(dtexNew$cpd.name, "exchange"))
  cmod@react_comp <- c(cmod@react_comp, rep("e",nrow(dtexNew)))
  cmod@react_attr <- rbind(cmod@react_attr,
                           data.frame(CVTerms = rep("",nrow(dtexNew)),
                                      SBOTerm = rep("http://identifiers.org/SBO:0000627",
                                                    nrow(dtexNew))))

  cmod@lowbnd <- c(cmod@lowbnd, dtexNew$lb)
  cmod@uppbnd <- c(cmod@uppbnd, dtexNew$ub)

  cmod@met_id <- c(cmod@met_id, dtexNew$cpd)
  cmod@met_name <- c(cmod@met_name, dtexNew$cpd.name)
  cmod@met_comp <- c(cmod@met_comp, rep("e",nrow(dtexNew)))
  cmod@met_attr <- rbind(cmod@met_attr,
                         dtexNew[,colnames(models[[1]]@met_attr)])

  #
  # Attach organism-specific exchange reactions to shared community
  # extracellular space
  #
  orgExInds <- grep("^M[0-9]+_EX_",cmod@react_id)
  tmpExS <- cmod@S[,orgExInds]
  excpds <- gsub("^M[0-9]+_","",cmod@met_id[tmpExS@i+1])
  cpdExInds <- match(excpds, cmod@met_id)

  cmod@S[matrix(c(cpdExInds,orgExInds), ncol = 2)] <- 1

  #
  # Open bounds if wanted
  #
  if(origExchangeBounds %in% c("open","open_lb")) {
    cmod <- changeBounds(cmod, grep("^M[0-9]+_EX_",cmod@react_id),
                         lb = -COBRAR_SETTINGS("MAXIMUM"))
  }
  if(origExchangeBounds %in% c("open","open_ub")) {
    cmod <- changeBounds(cmod, grep("^M[0-9]+_EX_",cmod@react_id),
                         ub = COBRAR_SETTINGS("MAXIMUM"))
  }

  #
  # Genes and GPRs
  #
  cmod@allGenes <- unlist(lapply(indeces,
                                 function(idx) paste0(idx,"_",
                                                      models[[idx]]@allGenes)))
  genetmp <- lapply(models, function(x) x@genes)
  cmod@genes <- unlist(lapply(names(genetmp), function(x) lapply(genetmp[[x]], function(y, idx) {
    if(length(y) > 0)
      return(paste0(idx,"_",y))
  }, idx = x)), recursive = FALSE)
  cmod@gprRules <- unlist(lapply(models, function(x) x@gprRules))
  cmod@genes_attr <- do.call("rbind",lapply(models, function(x) x@genes_attr))

  #
  # Subsystems
  #

  # TODO

  #
  # Constraints
  #
  nconstr <- unlist(lapply(models, constraint_num))
  constroffset <- cumsum(c(0, nrxns[-n]))
  cmod@constraints@coeff <- Matrix(0,
                                   nrow = sum(nconstr),
                                   ncol = ncol(cmod@S),
                                   sparse = TRUE)

  inds <- lapply(models, function(x)  {
    tmp <- which(x@constraints@coeff!=0, arr.ind = TRUE)
    tmp <- cbind(tmp, x@constraints@coeff[tmp])
  })

  for(i in 1:length(inds)) {
    indstmp <- inds[[i]]
    indstmp[,"row"] <- indstmp[,"row"]+constroffset[i]
    indstmp[,"col"] <- indstmp[,"col"]+rxnoffset[i]
    cmod@constraints@coeff[indstmp[,c("row","col")]] <- indstmp[,3]
  }

  cmod@constraints@lb <- unlist(lapply(models, function(x) x@constraints@lb))
  cmod@constraints@ub <- unlist(lapply(models, function(x) x@constraints@ub))
  cmod@constraints@rtype <- unlist(lapply(models, function(x) x@constraints@rtype))

  #
  # Objective coefficients
  #
  cmod@obj_coef <- unlist(lapply(models, function(x) x@obj_coef))
  cmod@obj_coef <- c(cmod@obj_coef, rep(0, ncol(cmod@S)-length(cmod@obj_coef)))

  return(cmod)
}
