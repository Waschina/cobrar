#------------------------------------------------------------------------------#
# Functions to modify the structure / properties of a ModelOrg                 #
#------------------------------------------------------------------------------#

#' Change flux bounds
#'
#' The function changes either upper bounds, lower bounds, or both for specific
#' reactions.
#'
#' @param model Model of class \link{ModelOrg}
#' @param react A character vector stating the reaction IDs in a model or a
#' numeric vector providing the reaction indices.
#' @param lb A numeric vector giving the new lower flux bounds for reactions
#' \code{react}. If \code{lb} is of length 1, the same value will be used for
#' all reactions.
#' @param ub A numeric vector giving the new upper flux bounds for reactions
#' \code{react}. If \code{ub} is of length 1, the same value will be used for
#' all reactions.
#'
#' @returns An updated model of class \link{ModelOrg}
#'
#' @family Model manipulation tools
#' @export
changeBounds <- function(model, react, lb = NULL, ub = NULL) {

  stopifnot("lb must be numeric" = is.null(lb) || is.numeric(lb),
            "ub must be numeric" = is.null(ub) || is.numeric(ub))

  if(!all(checkReactId(model, react))) {
    stop("Please check your reaction IDs/indices in argument 'react'.")
  }

  react.idx <- react_pos(model, react)

  # the actual change
  if(!is.null(lb)) {
    model@lowbnd[react.idx] <- lb
  }
  if(!is.null(ub)) {
    model@uppbnd[react.idx] <- ub
  }

  return(model)
}

#' Remove reactions from a model
#'
#' This function removes specified reactions from a model.
#'
#' @param model Model of class \link{ModelOrg}
#' @param react A character vector stating the reaction IDs in a model or a
#' numeric vector providing the reaction indices.
#' @param rm_met Logical. Should metabolites, which are singletons after the
#' reaction removal, be deleted as well?
#'
#' @returns An updated model of class \link{ModelOrg}
#'
#' @note
#' If the reaction participates in a user constraint, this constraint is
#' removed from the model.
#'
#' @family Model manipulation tools
#' @export
rmReact <- function(model, react, rm_met = TRUE) {
  if(length(react) == 0)
    return(model)

  if(!all(checkReactId(model, react))) {
    stop("Please check your reaction IDs/indices in argument 'react'.")
  }

  react <- react_pos(model, react)

  model@S        <- model@S[,-react, drop = FALSE]
  model@obj_coef <- model@obj_coef[-react]
  model@subSys   <- model@subSys[-react,, drop = FALSE]

  model@react_id   <- model@react_id[-react]
  model@react_name <- model@react_name[-react]
  model@react_comp <- model@react_comp[-react]
  model@lowbnd     <- model@lowbnd[-react]
  model@uppbnd     <- model@uppbnd[-react]
  model@react_attr <- model@react_attr[-react,, drop = FALSE]

  model@gprRules <- model@gprRules[-react]
  model@genes    <- model@genes[-react]

  rmconstr <- which(model@constraints@coeff[,react, drop = FALSE] != 0, arr.ind = T)[,1]
  if(length(rmconstr) > 0)
    model <- rmConstraint(model, rmconstr)

  model@constraints@coeff <- model@constraints@coeff[,-react, drop = FALSE]

  if(rm_met) {
    metrm <- which(apply(model@S,1,function(x) all(x == 0))) # identifies unused mets
    if(length(metrm) > 0)
      model <- rmMetabolite(model, metrm)
  }

  return(model)
}

#' Remove genes from a model
#'
#' This function removes specified genes from a model, and optionally also
#' reactions and metabolites inaccessible after gene knock outs.
#'
#' @param model Model of class \link{ModelOrg}
#' @param gene A character vector stating the reaction IDs in a model or a
#' numeric vector providing the reaction indices.
#' @param rm_react Logical. Should reaction, which are inaccessible after the
#' gene knock outs, be deleted as well?
#' @param rm_met Logical. Should metabolites, which are singletons after the
#' reaction removal, be deleted as well?
#'
#' @returns An updated model of class \link{ModelOrg}
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#' mod
#'
#' # create a double gene knock-out mutant
#' mod_KO <- rmGene(mod, c("b4152","b0116"))
#' mod_KO
#'
#' @family Model manipulation tools
#' @export
rmGene <- function(model, gene, rm_react = TRUE, rm_met = TRUE) {
  if(length(gene) == 0)
    return(model)

  if(!all(checkGeneId(model, gene))) {
    stop("Please check your gene IDs/indices in argument 'gene'.")
  }

  gene_ids <- gene
  gene <- gene_pos(model, gene)

  rmReactions <- geneDel(model, gene)

  # rm gene parts
  model@allGenes <- model@allGenes[-gene]
  model@genes_attr <- model@genes_attr[-gene,]
  model@genes <- lapply(model@genes,
                        function(x) ifelse(x %in% gene_ids, NA_character_,x))

  # rm reaction (and metabolites)
  if(rm_react)
    model <- rmReact(model, rmReactions, rm_met)

  return(model)
}

#' Add constraints to model
#'
#' Add linear reaction flux constraints to a metabolic network.
#'
#' @param model Model of class \link{ModelOrg}
#' @param react Character vector or a list of character vectors containing the
#' model's reactions IDs that are part of the respective constraint.
#' @param coeff Numeric vector or list of numeric vectors defining the
#' coefficients for the reactions listed in 'react'.
#' @param rtype Character vector describing the type of the constraint(s). See
#' details.
#' @param lb,ub Numeric vector defining the lower and upper bound(s) of the
#' constraint(s).
#'
#' @details
#' The slot "rtype" describes the type of each constraint. Valid values and
#' their effects are:
#' | *code* | *description* | *rule* |
#' | :----: | :--- | :----: |
#' | "F" | free constraint | \eqn{-\infty < x < \infty} |
#' | "L" | constraint with lower bound | \eqn{lb \leq x \leq \infty} |
#' | "U" | constraint with upper bound | \eqn{-\infty \leq x \leq ub} |
#' | "D" | double-bounded (ranged) constraint | \eqn{lb \leq x \leq ub} |
#' | "E" | fixed (equality constraint) | \eqn{lb = x = ub} |
#'
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#'
#' # Simulate anaerobic growth
#' mod <- changeBounds(mod, "EX_o2_e", lb = 0)
#'
#' # Limit the proton production depending on the growth rate
#' mod <- addConstraint(mod,
#'                      react = c("EX_h_e","BIOMASS_Ecoli_core_w_GAM"),
#'                      coeff = c(1, -20), rtype = "U", ub = 0)
#'
#'
#' @importFrom stats var
#'
#' @docType methods
#' @rdname addConstraint-methods
#'
#' @family Model manipulation tools
#' @export
setGeneric("addConstraint" ,valueClass = "ModelOrg", function(model,
                                                              react,
                                                              coeff,
                                                              rtype,
                                                              lb = NULL,
                                                              ub = NULL) {
  standardGeneric("addConstraint")
})
#' @rdname addConstraint-methods
#' @aliases addConstraint,ModelOrg,character,numeric,character
setMethod("addConstraint", signature(model = "ModelOrg",
                                     react = "character",
                                     coeff = "numeric",
                                     rtype = "character"),
          function(model, react, coeff, rtype, lb = NULL, ub = NULL) {
            return(addConstraint(model,
                                 react = list(react),
                                 coeff = list(coeff),
                                 lb = lb,
                                 ub = ub,
                                 rtype = rtype))
          }
)
#' @rdname addConstraint-methods
#' @aliases addConstraint,ModelOrg,list,list,character
setMethod("addConstraint", signature(model = "ModelOrg",
                                     react = "list",
                                     coeff = "list",
                                     rtype = "character"),
          function(model, react, coeff, rtype, lb = NULL, ub = NULL) {

            nc <- length(react)

            if(is.null(lb))
              lb <- rep(NA_real_, nc)
            if(is.null(ub))
              ub <- rep(NA_real_, nc)

            # validity checks
            if(var(c(length(react),
                     length(coeff),
                     length(rtype),
                     length(lb),
                     length(ub))) != 0) {
              "Lengths all arguments 'react', 'coeff', 'rtype', 'lb', and 'ub' must be equal."
            }
            if(any(!(rtype %in% c("F","L","U","D","E"))))
              stop("Elements of the vector 'rtype' must be \"F\", \"L\", \"U\", \"D\", or \"E\".")
            if(any(rtype != "F" & is.na(lb) & is.na(ub)))
              stop("'ub' and 'lb' can only be both undefined if 'rtype' is \"F\".")
            if(any(rtype == "L" & is.na(lb)))
              stop("If 'rtype' is \"L\", 'lb' cannot be undefined.")
            if(any(rtype == "U" & is.na(ub)))
              stop("If 'rtype' is \"U\", 'ub' cannot be undefined.")
            if(any(rtype == "D" & (is.na(ub) | is.na(lb))))
              stop("If 'rtype' is \"D\", 'ub' and 'lb' both need to be defined.")

            indtmp <- which(rtype == "E" & is.na(lb))
            lb[indtmp] <- ub[indtmp]
            indtmp <- which(rtype == "E" & is.na(ub))
            ub[indtmp] <- lb[indtmp]

            if(any(unlist(lapply(react, length)) != unlist(lapply(coeff, length))))
              stop("List elementes of 'react' must have the same length as the corresponding list elements in 'coeff'.")

            if(any(unlist(lapply(react, duplicated))))
              stop("'react' IDs cannot be duplicated within a constraint definition.")

            if(any(unlist(lapply(react, function(x) !(x %in% model@react_id)))))
              stop("Not all reaction IDs in 'react' are part of the model.")

            I <- matrix(c(rep(1:nc, unlist(lapply(react, length))),
                          unlist(lapply(react, function(x) match(x, model@react_id)))),
                        ncol = 2)


            out <- Matrix(0, nrow = nc, ncol = react_num(model), sparse = T)
            out[I] <- unlist(coeff)

            model@constraints@coeff <- rbind(model@constraints@coeff,
                                              out)
            model@constraints@lb <- c(model@constraints@lb, lb)
            model@constraints@ub <- c(model@constraints@ub, ub)
            model@constraints@rtype <- c(model@constraints@rtype, rtype)

            return(rmDuplicateConstraints(model))
          }
)

#' Remove constraints
#'
#' Remove specific user constraints from a metabolic model.
#'
#' @param model Model of class \link{ModelOrg}
#' @param ind Integer vector with the indices of the constraints to be removed
#'
#' @seealso [printConstraint()]
#'
#' @family Model manipulation tools
#' @export
rmConstraint <- function(model, ind) {
  if(length(ind) == 0)
    return(model)

  if(constraint_num(model) == 0 || any(!(ind %in% 1:constraint_num(model)))) {
    stop("Invalid index or indices for constraints.")
  }

  ind <- unique(ind)

  model@constraints@coeff <- model@constraints@coeff[-ind,, drop = FALSE]
  model@constraints@lb <- model@constraints@lb[-ind]
  model@constraints@ub <- model@constraints@ub[-ind]
  model@constraints@rtype <- model@constraints@rtype[-ind]

  return(model)
}

#' Add or modify a reaction
#'
#' The function can be used to add or modify a reaction in an existing model.
#'
#' @param model Model of class \link{ModelOrg}
#' @param id Character for the reaction ID
#' @param met Character vector providing the IDs of metabolites that participate
#' in the reaction
#' @param Scoef Numeric vector (same length as `met`) of stoichiometric
#' coefficients for the metabolites in `met`. The value in `Scoef[i]` is the
#' stoichiometric coefficient of the metabolite in `met[i]`.
#' @param reversible This option has now effect and is only here for legacy
#' reasons. Whether a reaction is reversible or not is inferred by cobrar based
#' on the lower and upper bounds.
#' @param lb,ub Single numeric values that define the lower and upper flux
#' limits, respectively.
#' @param obj Single numeric value for the coefficient of the reaction in the
#' objective function.
#' @param subsystem A vector of character strings containing the sub system IDs
#' to which the reaction belongs.
#' @param subsystemName A character vector (same length as `subsystem`) for the
#' names of the subsystems. If the subsystem is already part of the model and
#' you do not want to change its name, just use NA the corresponding entry.
#' @param gprAssoc A single character string giving the Gene-Product-Reaction
#' (GPR) association for the reaction. If NA, no GRP association is created.
#' @param reactName A single character string giving the name for the reaction.
#' If NA, the value of argument `id` is used.
#' @param metName A vector of character strings of the same length as `met`
#' containing the metabolites names for the metabolites given in argument `met`.
#' @param metComp A vector of character strings of the same length as `met`
#' specifying the compartment IDs for the metabolites given in argument `met`.
#' @param metCharge A numeric vector of the same length as `met` defining the
#' charges for the metabolites given in argument `met`.
#' @param metChemicalFormula A character vector of the same length as `met`
#' defining the chemical formulas for the metabolites given in argument `met`.
#' @param CVTerms Cross-references to other resources.
#' @param SBOTerm A termID from the Systems Biology Ontology.
#'
#' @details
#' If you want to use the function to update data of a pre-existing reaction but
#' not its stoichiometry, use NA for the parameters 'met' and 'Scoeff'.
#' If the reaction is already part of the model, any reaction value (e.g., lb,
#' ub, reactName), that is set to NA has the effect that the old value will be
#' used.
#' If the reaction is already part of the model, and values for the parameter
#' `subsystem` are provided, all previous set subsystem associations of the
#' reaction will be removed.
#' If metabolites or subsystems are not part of the model yet, they will be
#' added.
#'
#' @examples
#' # This example adds the 4-aminobutyrate degradation pathway to the E. coli
#' # core metabolic model
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#'
#' fba(mod)
#'
#' # 4abut transport: 4abut_e + h_e <=> 4abut_c + h_c
#' mod <- addReact(mod, id = "ABUTt", Scoef = c(-1,-1,1,1),
#'                 met = c("4abut_e","h_e","4abut_c","h_c"), reversible = TRUE,
#'                 lb = -1000, ub = 1000,
#'                 reactName = "4-aminobutyrate transport in via proton symport",
#'                 metName = c("4-aminobutyrate",NA, "4-aminobutyrate",NA),
#'                 metComp = c("e","e","c","c"), metCharge = c(0,NA,0,NA),
#'                 metChemicalFormula = c("C4H9NO2",NA,"C4H9NO2",NA),
#'                 SBOTerm = "SBO:0000185")
#'
#' # exchange reaction for 4abut (with 1.5 mmol/gDW/hr availability)
#' mod <- addReact(mod, id = "EX_4abut_e", Scoef = c(-1), met = "4abut_e",
#'                 lb = -1.5, ub = 1000, reactName = "4-aminobutyrate exchange",
#'                 SBOTerm = "SBO:0000627")
#'
#' # 4abut amninotransferase (EC 2.6.1.19)
#' mod <- addReact(mod, id = "ABTA", Scoef = c(-1,-1,1,1),
#'                 met = c("4abut_c","akg_c","glu__L_c","sucsal_c"),
#'                 lb = 0,
#'                 reactName = "4-aminobutyrate transaminase",
#'                 metName = c(NA,NA,NA,"Succinic semialdehyde"),
#'                 metComp = c(NA,NA,NA,"c"), metCharge = c(NA,NA,NA,-1),
#'                 metChemicalFormula = c(NA,NA,NA,"C4H5O3"),
#'                 subsystem = "GABAdegr", subsystemName = "4-aminobutyrate degradation",
#'                 CVTerms = "bqbiol_is;http://identifiers.org/ec-code/2.6.1.19",
#'                 gprAssoc = "b2662 | b1302")
#'
#' # Succinate-semialdehyde dehydrogenase (NAD) (EC 1.2.1.24)
#' mod <- addReact(mod, id = "SSALx", Scoef = c(-1,-1,-1,2,1,1),
#'                 met = c("h2o_c","nad_c","sucsal_c","h_c","nadh_c","succ_c"),
#'                 lb = 0,
#'                 reactName = "Succinate-semialdehyde dehydrogenase (NAD)",
#'                 subsystem = "GABAdegr",
#'                 CVTerms = "bqbiol_is;http://identifiers.org/ec-code/1.2.1.24",
#'                 gprAssoc = "b1525")
#'
#' printReaction(mod, "SSALx")
#'
#' fba(mod)
#'
#' @family Model manipulation tools
#' @export
addReact <- function(model,
                     id,
                     met,
                     Scoef,
                     reversible = FALSE,
                     lb = 0,
                     ub = COBRAR_SETTINGS("MAXIMUM"),
                     obj = 0,
                     subsystem = NA,
                     subsystemName = NA,
                     gprAssoc = NA,
                     reactName = NA,
                     metName = NA,
                     metComp = NA,
                     metCharge = NA,
                     metChemicalFormula = NA,
                     CVTerms = NA,
                     SBOTerm = "SBO:0000176") {

  #--------------#
  # basic checks #
  #--------------#
  if(any(duplicated(met)))
    stop("Duplicates in metabolite IDs.")
  if(length(met) == 0)
    stop("Reaction needs at least one participating metabolite.")
  if(length(Scoef) != length(met))
    stop("Mismatch of number of metabolites and provides stoichiometrix coefficients.")

  #--------------------------------------#
  # Check if reaction addition or update #
  #--------------------------------------#
  if(id %in% model@react_id) {
    indR <- which(id == model@react_id)

    # remove previous stoichiometry if new coefficients are provided for
    # existing reaction:
    if(!is.na(Scoef[1]))
      model@S[,indR] <- rep(0, met_num(model))

  } else {
    indR <- react_num(model) + 1
    if(is.na(reactName))
      reactName <- id

    # extend data structures

    # S and obj.-coeff.
    newcol <- as(as(as(Matrix(data=0, nrow = nrow(model@S), ncol = 1,
                              sparse = TRUE),
                       "dMatrix"),
                    "generalMatrix"),
                 "CsparseMatrix")
    model@S <- cbind(model@S, newcol)
    model@obj_coef <- append(model@obj_coef, 0)

    # reaction slots
    model@react_attr[react_num(model)+1,] <- NA
    model@react_comp <- append(model@react_comp, NA)
    model@react_id   <- append(model@react_id, id)
    model@react_name <- append(model@react_name, NA_character_)

    # bounds
    model@lowbnd <- append(model@lowbnd, NA)
    model@uppbnd <- append(model@uppbnd, NA)

    # GPRs
    model@genes <- append(model@genes, list(character(0L)))
    model@gprRules <- append(model@gprRules, "")

    # constraints
    model@constraints@coeff <- cbind(model@constraints@coeff,
                                     matrix(0, nrow = constraint_num(model),
                                            ncol = 1))

    # subsystems
    model@subSys <- rbind(model@subSys,
                          Matrix(FALSE, nrow = 1, ncol = ncol(model@subSys),
                                 sparse = TRUE))
  }

  # Add/update metabolites if necessary
  model <- addMetabolite(model = model, id = met, name = metName, comp = metComp,
                         chemicalFormula = metChemicalFormula, charge = metCharge)
  indMs <- match(met, model@met_id)

  # Add/update subsystems  if necessary
  if(length(subsystem) > 0 && !any(is.na(subsystem))) {
    model <- addSubsystem(model, id = subsystem, name = subsystemName)
    indSubSys <- subsys_pos(model, subsystem)
  }

  #----------------------------#
  # Add/update reaction values #
  #----------------------------#

  # S and obj.-coeff.
  if(!is.na(Scoef[1]))
    model@S[matrix(c(indMs, rep(indR,length(indMs))),ncol = 2)] <- Scoef
  if(!is.na(obj))
    model@obj_coef[indR] <- obj

  # reaction slots
  if(!is.na(reactName))
    model@react_name[indR] <- reactName
  if(!is.na(CVTerms))
    model@react_attr$CVTerms[indR] <- CVTerms
  if(!is.na(SBOTerm))
    model@react_attr$SBOTerm[indR] <- SBOTerm

  # bounds
  model@lowbnd[indR] <- lb
  model@uppbnd[indR] <- ub

  # subsys
  if(length(subsystem) > 0 && !any(is.na(subsystem))) {
    model@subSys[indR,] <- FALSE
    model@subSys[matrix(c(rep(indR,length(indSubSys)), indSubSys),ncol = 2)] <- TRUE
  }

  # GPR
  if(!is.na(gprAssoc)) {
    gpr_new <- parseBoolean(gprAssoc)
    model@genes[indR] <- list(gpr_new$gene)
    model@gprRules[indR] <- gpr_new$rule

    # in case the GPR involves new genes
    newgenes <- gpr_new$gene[!(gpr_new$gene %in% model@allGenes) & gpr_new$gene != ""]
    if(length(newgenes) > 0)
      model <- addGene(model, newgenes)
  }

  return(model)
}

#' Add metabolites or update their data
#'
#' The function allows you to add one or more metabolites to a model. When
#' providing the ID of an already existing metabolite, you can use this function
#' to update metabolite information.
#'
#' @param model Model of class \link{ModelOrg}
#' @param id Character vector with metabolite IDs
#' @param comp Character vector of the metabolites' compartment IDs
#' @param name Character vector for metabolite names
#' @param chemicalFormula Character vector for the metabolites' chemical formulas
#' @param charge Numeric vector for the metabolites' charge
#' @param CVTerms Character vector for the metabolites' CV-Terms
#' @param SBOTerm Character vector for the metabolites' SBO-Terms
#'
#' @family Model manipulation tools
#' @export
addMetabolite <- function(model, id, comp = NA, name = NA, chemicalFormula = NA,
                          charge = NA, CVTerms = NA,
                          SBOTerm = rep("SBO:0000247", length(id))) {
  norig <- met_num(model)

  #--------------------------------#
  # If optional values are missing #
  #--------------------------------#
  if(length(comp) == 1 && is.na(comp))
    comp <- rep(NA, length(id))
  if(length(name) == 1 && is.na(name))
    name <- rep(NA, length(id))
  if(length(chemicalFormula) == 1 && is.na(chemicalFormula))
    chemicalFormula <- rep(NA, length(id))
  if(length(charge) == 1 && is.na(charge))
    charge <- rep(NA, length(id))
  if(length(CVTerms) == 1 && is.na(CVTerms))
    CVTerms <- rep(NA, length(id))
  if(length(SBOTerm) == 1 && is.na(SBOTerm))
    SBOTerm <- rep(NA, length(id))

  #--------------#
  # basic checks #
  #--------------#
  if(any(duplicated(id)))
    stop("Duplicates in metabolite IDs.")
  if(length(id) == 0)
    stop("No metabolite ID provided.")
  if(length(name) != length(id))
    stop("Mismatch of number of metabolite IDs and Names.")
  if(length(comp) != length(id))
    stop("Mismatch of number of metabolite IDs and compartment.")
  if(length(chemicalFormula) != length(id))
    stop("Mismatch of number of metabolite IDs and chemical formulas.")
  if(length(charge) != length(id))
    stop("Mismatch of number of metabolite IDs and charge.")


  #---------------------------------------------------------------#
  # check if metabolites are new or if updated infos are provided #
  #---------------------------------------------------------------#
  indM <- match(id, model@met_id)
  nnew <- sum(is.na(indM))
  if(nnew > 0) {
    indnew <- which(is.na(indM))
    indM[is.na(indM)] <- (norig+1):(norig+nnew)

    # extent metabolite-related model structures

    # S
    newrows <- as(as(as(Matrix(data=0, nrow = nnew, ncol = ncol(model@S), sparse = TRUE),
                        "dMatrix"),
                     "generalMatrix"),
                  "CsparseMatrix")
    model@S <- rbind(model@S, newrows)

    # metabolite slots
    model@met_attr[norig+nnew,] <- NA
    model@met_comp <- append(model@met_comp, rep(NA, nnew))
    model@met_name <- append(model@met_name, rep(NA_character_, nnew))
    model@met_id   <- append(model@met_id, id[indnew])

    # use ids if names for new metabolites are not provided
    name[which(indM > norig & is.na(name))] <- id[which(indM > norig & is.na(name))]

  }

  #------------------------------#
  # Add/update metabolite values #
  #------------------------------#

  # name
  model@met_name[indM][which(!is.na(name))] <- name[which(!is.na(name))]

  # attributes
  model@met_attr$chemicalFormula[indM][which(!is.na(chemicalFormula))] <- chemicalFormula[which(!is.na(chemicalFormula))]
  model@met_attr$charge[indM][which(!is.na(charge))] <- charge[which(!is.na(charge))]
  model@met_attr$CVTerms[indM][which(!is.na(CVTerms))] <- CVTerms[which(!is.na(CVTerms))]
  model@met_attr$SBOTerm[indM][which(!is.na(SBOTerm))] <- SBOTerm[which(!is.na(SBOTerm))]

  # compartment
  comp_new <- comp_pos(model, comp)
  if(any(is.na(comp_new[!is.na(comp)]))) {
    model <- addCompartment(model, unique(comp[is.na(comp_new) & !is.na(comp)]))
    comp_new <- comp_pos(model, comp)
  }
  model@met_comp[indM][which(!is.na(comp))] <- model@mod_compart[comp_new[which(!is.na(comp))]]


  return(model)
}

#' Remove metabolites from a model
#'
#' This function removes specified metabolites from a model.
#'
#' @param model Model of class \link{ModelOrg}
#' @param met A character vector stating the metabolite IDs in a model or a
#' numeric vector providing the metabolite indices.
#'
#' @returns An updated model of class \link{ModelOrg}
#'
#' @note
#' If at least one of the provided metabolites still participates in a reaction,
#' the function stops with an error message.
#'
#' @family Model manipulation tools
#' @export
rmMetabolite <- function(model, met) {
  if(length(met) == 0)
    return(model)

  if(!all(checkMetId(model, met))) {
    stop("Please check your metabolite IDs/indices in argument 'met'.")
  }

  met <- met_pos(model, met)

  # print(paste0("Removing metabolites: ",paste(model@met_id[met], collapse = ", ")))

  if(any(model@S[met,] != 0)) {
    stop("At least one provided metabolite still participates in at least one reaction.")
  }

  # remove metabolite from data structures
  model@S        <- model@S[-met,, drop = FALSE]
  model@met_id   <- model@met_id[-met]
  model@met_name <- model@met_name[-met]
  model@met_comp <- model@met_comp[-met]
  model@met_attr <- model@met_attr[-met,, drop = FALSE]

  return(model)
}



#' Add genes or update their data
#'
#' The function allows you to add one or more genes to a model. When
#' providing the ID of an already existing genes, you can use this function
#' to update the gene's information.
#'
#' @param model Model of class \link{ModelOrg}
#' @param id Character vector with gene IDs
#' @param name Character vector for gene names
#' @param CVTerms Character vector for the genes' CV-Terms
#' @param SBOTerm Character vector for the genes' SBO-Terms
#'
#' @family Model manipulation tools
#' @export
addGene <- function(model, id, name = NA, CVTerms = NA,
                    SBOTerm = rep("SBO:0000243",length(id))) {
  norig <- gene_num(model)

  #--------------------------------#
  # If optional values are missing #
  #--------------------------------#
  if(length(name) == 1 && is.na(name))
    name <- rep(NA, length(id))
  if(length(CVTerms) == 1 && is.na(CVTerms))
    CVTerms <- rep(NA, length(id))
  if(length(SBOTerm) == 1 && is.na(SBOTerm))
    SBOTerm <- rep(NA, length(id))

  #--------------#
  # basic checks #
  #--------------#
  if(any(duplicated(id)))
    stop("Duplicates in gene IDs.")
  if(length(id) == 0)
    stop("No gene ID provided.")
  if(length(name) != length(id))
    stop("Mismatch of number of gene IDs and Names.")
  if(length(CVTerms) != length(id))
    stop("Mismatch of number of gene IDs and CVTerms.")
  if(length(SBOTerm) != length(id))
    stop("Mismatch of number of gene IDs and SBOTerms.")

  #---------------------------------------------------------#
  # check if genes are new or if updated infos are provided #
  #---------------------------------------------------------#
  indG <- match(id, model@allGenes)
  nnew <- sum(is.na(indG))
  if(nnew > 0) {
    indnew <- which(is.na(indG))
    indG[is.na(indG)] <- (norig+1):(norig+nnew)

    # extent gene-related model structures

    # gene slots
    model@genes_attr[norig+nnew,] <- NA
    model@allGenes   <- append(model@allGenes, id[indnew])

    # use ids if names for new genes are not provided
    name[which(indG > norig & is.na(name))] <- id[which(indG > norig & is.na(name))]

  }

  #------------------------#
  # Add/update gene values #
  #------------------------#

  # attributes
  model@genes_attr$name[indG][which(!is.na(name))] <- name[which(!is.na(name))]
  model@genes_attr$CVTerms[indG][which(!is.na(CVTerms))] <- CVTerms[which(!is.na(CVTerms))]
  model@genes_attr$SBOTerm[indG][which(!is.na(SBOTerm))] <- SBOTerm[which(!is.na(SBOTerm))]


  return(model)
}

#' Add compartments or update their data
#'
#' The function allows you to add one or more compartments to a model. When
#' providing the ID of an already existing compartment, you can use this
#' function to update the compartment's name.
#'
#' @param model Model of class \link{ModelOrg}
#' @param id Character vector with compartment IDs
#' @param name Character vector for compartment names
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#' mod <- addCompartment(mod, id = "p", name = "periplasm")
#'
#' @family Model manipulation tools
#' @export
addCompartment <- function(model, id, name = NA) {
  norig <- comp_num(model)

  #--------------------------------#
  # If optional values are missing #
  #--------------------------------#
  if(length(name) == 1 && is.na(name))
    name <- rep(NA, length(id))

  #--------------#
  # basic checks #
  #--------------#
  if(any(duplicated(id)))
    stop("Duplicates in compartment IDs.")
  if(length(id) == 0)
    stop("No compartment ID provided.")
  if(length(name) != length(id))
    stop("Mismatch of number of compartment IDs and Names.")

  #----------------------------------------------------------------#
  # check if compartments are new or if updated infos are provided #
  #----------------------------------------------------------------#
  indC <- comp_pos(model, id)
  nnew <- sum(is.na(indC))
  if(nnew > 0) {
    indnew <- which(is.na(indC))
    indC[is.na(indC)] <- (norig+1):(norig+nnew)

    # extent compartment-related model structures
    model@mod_compart      <- append(model@mod_compart, id[indnew])
    model@mod_compart_name <- append(model@mod_compart_name, name[indnew])

    # use ids if names for new metabolites are not provided
    name[which(indC > norig & is.na(name))] <- id[which(indC > norig & is.na(name))]
  }

  #-------------------------------#
  # Add/update compartment values #
  #-------------------------------#

  # name
  model@mod_compart_name[indC][which(!is.na(name))] <- name[which(!is.na(name))]

  return(model)
}

#' Remove compartments from a model
#'
#' This function removes specified compartments from a model.
#'
#' @param model Model of class \link{ModelOrg}
#' @param comp A character vector stating the compartment IDs in a model or a
#' numeric vector providing the compartment indices.
#'
#' @returns An updated model of class \link{ModelOrg}
#'
#' @note
#' If at least one of the provided compartments still has metabolites associated
#' with it, the function stops with an error message.
#'
#' @family Model manipulation tools
#' @export
rmCompartment <- function(model, comp) {
  if(length(comp) == 0)
    return(model)

  if(!all(checkCompartmentId(model, comp))) {
    stop("Please check your metabolite IDs/indices in argument 'comp'.")
  }

  comp <- comp_pos(model, comp)

  if(any(comp_pos(model, model@met_comp) %in% comp)) {
    stop("At least one provided compartment still has metabolites associated with it.")
  }

  # remove metabolite from data structures
  model@mod_compart      <- model@mod_compart[-comp]
  model@mod_compart_name <- model@mod_compart_name[-comp]

  return(model)
}

#' Add subsystems or update their data
#'
#' The function allows you to add one or more subsystems to a model. When
#' providing the ID of an already existing subsystem, you can use this
#' function to update the subsystem's name.
#'
#' @param model Model of class \link{ModelOrg}
#' @param id Character vector with subsystem IDs
#' @param name Character vector for subsystem names
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#' mod <- addSubsystem(mod, id = "Bifidoshunt",
#'                     name = "glucose fermentation to acetate and lactate (Bifidobacteria)")
#'
#' @family Model manipulation tools
#' @export
addSubsystem <- function(model, id, name = NA) {
  norig <- subsys_num(model)

  #--------------------------------#
  # If optional values are missing #
  #--------------------------------#
  if(length(name) == 1 && is.na(name))
    name <- rep(NA, length(id))

  #--------------#
  # basic checks #
  #--------------#
  if(any(duplicated(id)))
    stop("Duplicates in subsystem IDs.")
  if(length(id) == 0)
    stop("No subsystem ID provided.")
  if(length(name) != length(id))
    stop("Mismatch of number of subsystem IDs and Names.")

  #----------------------------------------------------------------#
  # check if subsystems are new or if updated names are provided   #
  #----------------------------------------------------------------#
  indS <- subsys_pos(model, id)
  nnew <- sum(is.na(indS))
  if(nnew > 0) {
    indnew <- which(is.na(indS))
    indS[is.na(indS)] <- (norig+1):(norig+nnew)

    # extent compartment-related model structures
    model@subSys_id   <- append(model@subSys_id, id[indnew])
    model@subSys_name <- append(model@subSys_name, rep(NA_character_,nnew))
    model@subSys <- cbind(model@subSys,
                          Matrix(FALSE,
                                 nrow = nrow(model@subSys),
                                 ncol = length(indnew),
                                 sparse = TRUE))
    colnames(model@subSys)[indnew] <- id[indnew]

    # use ids if names for new metabolites are not provided
    name[which(indS > norig & is.na(name))] <- id[which(indS > norig & is.na(name))]
  }

  #-------------------------------#
  # Add/update compartment values #
  #-------------------------------#

  # name
  model@subSys_name[indS][which(!is.na(name))] <- name[which(!is.na(name))]

  return(model)
}

#' Remove subsystems from a model
#'
#' This function removes specified subsystems from a model.
#'
#' @param model Model of class \link{ModelOrg}
#' @param subsystem A character vector stating the subsystem IDs in a model or a
#' numeric vector providing the subsystem indices.
#'
#' @returns An updated model  of class \link{ModelOrg}
#'
#' @family Model manipulation tools
#' @export
rmSubsystem <- function(model, subsystem) {
  if(length(subsystem) == 0)
    return(model)

  if(!all(checkSubsystemId(model, subsystem))) {
    stop("Please check your subsystem IDs/indices in argument 'subsystem'.")
  }

  subsystem <- subsys_pos(model, subsystem)

  # remove metabolite from data structures
  model@subSys      <- model@subSys[,-subsystem]
  model@subSys_id   <- model@subSys_id[-subsystem]
  model@subSys_name <- model@subSys_name[-subsystem]

  return(model)
}
