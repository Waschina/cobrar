#------------------------------------------------------------------------------#
# Functions to modify the structure / properties of a modelorg                 #
#------------------------------------------------------------------------------#

#' Change flux bounds
#'
#' The function changes either upper bounds, lower bounds, or both for specific
#' reactions.
#'
#' @param model Model of class \link{modelorg}
#' @param react A character vector stating the reaction IDs in a model or a
#' numeric vector providing the reaction indices.
#' @param lb A numeric vector giving the new lower flux bounds for reactions
#' \code{react}. If \code{lb} is of length 1, the same value will be used for
#' all reactions.
#' @param ub A numeric vector giving the new upper flux bounds for reactions
#' \code{react}. If \code{ub} is of length 1, the same value will be used for
#' all reactions.
#'
#' @returns An updated model of class \link{modelorg}
#'
#' @export
changeBounds <- function(model, react, lb = NULL, ub = NULL) {

  stopifnot("lb must be numeric" = is.null(lb) || is.numeric(lb),
            "ub must be numeric" = is.null(ub) || is.numeric(ub))

  if(!all(checkReactId(model, react))) {
    stop("Please check your reaction IDs/indices in argument 'react'.")
  }

  # the actual change
  if(!is.null(lb)) {
    model@lowbnd[match(react, model@react_id)] <- lb
  }
  if(!is.null(ub)) {
    model@lowbnd[match(react, model@react_id)] <- ub
  }

  return(model)
}

#' Remove reactions from a model
#'
#' This function removes specified reactions from a model.
#'
#' @param model Model of class \link{modelorg}
#' @param react A character vector stating the reaction IDs in a model or a
#' numeric vector providing the reaction indices.
#' @param rm_met Logical. Should metabolites, which are singletons after the
#' reaction removal, be deleted as well?
#'
#' @returns An updated model of class \link{modelorg}
#'
#' @note
#' If the reaction participates in a user constraint, this constraint is
#' removed from the model.
#'
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
  print(react)
  rmconstr <- which(model@constraints@coeff[,react, drop = FALSE] != 0, arr.ind = T)[,1]
  print(rmconstr)
  model <- rmConstraint(model, rmconstr)

  if(rm_met) {
    TMPmat <- as(model@S[,-react], "TsparseMatrix")
    metrm <- which(!(1:met_num(model) %in% (TMPmat@i+1)))
    if(length(metrm) > 0) {

      model@S        <- model@S[-metrm,, drop = FALSE]
      model@met_id   <- model@met_id[-metrm]
      model@met_name <- model@met_name[-metrm]
      model@met_comp <- model@met_comp[-metrm]
      model@met_attr <- model@met_attr[-metrm,, drop = FALSE]

    }
  }

  return(model)
}

#' Remove genes from a model
#'
#' This function removes specified genes from a model, and optionally also
#' reactions and metabolites inaccessible after gene knock outs.
#'
#' @param model Model of class \link{modelorg}
#' @param gene A character vector stating the reaction IDs in a model or a
#' numeric vector providing the reaction indices.
#' @param rm_react Logical. Should reaction, which are inaccessible after the
#' gene knock outs, be deleted as well?
#' @param rm_met Logical. Should metabolites, which are singletons after the
#' reaction removal, be deleted as well?
#'
#' @returns An updated model of class \link{modelorg}
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
#'
#' @export
rmGene <- function(model, gene, rm_react = TRUE, rm_met = TRUE) {
  if(length(gene) == 0)
    return(model)

  if(!all(checkGeneId(model, gene))) {
    stop("Please check your gene IDs/indices in argument 'gene'.")
  }

  gene <- gene_pos(model, gene)

  rmReactions <- geneDel(model, gene)

  # rm gene parts
  model@allGenes <- model@allGenes[-gene]
  model@allGenes_name <- model@allGenes_name[-gene]
  model@genes <- lapply(model@genes,
                        function(x) ifelse(x %in% gene, NA_character_,x))

  # rm reaction (and metabolites)
  if(rm_react)
    model <- rmReact(model, rmReactions)

  return(model)
}

#' Add constraints to model
#'
#' Add linear reaction flux constraints to a metabolic network.
#'
#' @param model Model of class \link{modelorg}
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
#' @export
setGeneric("addConstraint" ,valueClass = "modelorg", function(model,
                                                              react,
                                                              coeff,
                                                              rtype,
                                                              lb = NULL,
                                                              ub = NULL) {
  standardGeneric("addConstraint")
})
setMethod("addConstraint", signature(model = "modelorg",
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
setMethod("addConstraint", signature(model = "modelorg",
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
#' @param model Model of class \link{modelorg}
#' @param ind Integer vector with the indices of the constraints to be removed
#'
#' @seealso [printConstraint()]
#'
#' @export
rmConstraint <- function(model, ind) {
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
#' @param model Model of class \link{modelorg}
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
#' @param subSystem A vector of character strings containing the sub system IDs
#' to which the reaction belongs.
#' @param subSystemName A character vector (same length as `subSystem`) for the
#' names of the subsystems. If the subsystem is already part of the model and
#' you do not want to change its name, just use NA the corresponding entry.
#' @param gprAssoc A single character string giving the Gene-Product-Reaction
#' (GPR) association for the reaction. If NA, no GRP association is created.
#' @param reactName A single character string giving the name for the reaction.
#' If NA, the value of argument `id` is used.
#' @param metName A vector of character strings of the same length as `met`
#' containing the metabolites names for the metabolites given in argument `met`.
#' @param metCharge A numeric vector of the same length as `met` defining the
#' charges for the metabolites given in argument `met`.
#' @param metChemicalFormula A character vector of the same length as `met`
#' defining the chemical formulas for the metabolites given in argument `met`.
#' @param annotation An annotation string for the reaction.
#'
#' @details
#' If you want to use the function to update data of a pre-existing reaction but
#' not its stoichiometry, use NA for the parameters 'met' and 'Scoeff'.
#' If the reaction is already part of the model, any reaction value (e.g., lb,
#' ub, reactName), that is set to NA has the effect that the old value will be
#' used.
#' If the reaction is already part of the model, and values for the parameter
#' `subSystem` are provided, all previous set Subsystem associations of the
#' reaction will be removed.
#' If metabolites or subsystems are not part of the model yet, they will be
#' added.

#'
addReact <- function(model,
                     id,
                     met,
                     Scoef,
                     reversible = FALSE,
                     lb = 0,
                     ub = COBRAR_SETTINGS("MAXIMUM"),
                     obj = 0,
                     subSystem = NA,
                     subSystemName = NA,
                     gprAssoc = NA,
                     reactName = NA,
                     metName = NA,
                     metComp = NA,
                     metCharge = NA,
                     metChemicalFormula = NA,
                     annotation = NA) {

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
    model@S <- cbind(model@S, rep(0, met_num(model)))
    model@obj_coef <- append(model@obj_coef, 0)

    # reaction slots
    model@react_attr <- rbind(model@react_attr,
                            data.frame(annotation = NA_character_))
    model@react_comp <- append(model@react_comp, NA)
    model@react_id   <- append(model@react_id, id)
    model@react_name <- append(model@react_name, NA_character_)

    # bounds
    model@lowbnd <- append(model@lowbnd, NA)
    model@uppbnd <- append(model@uppbnd, NA)

    # GPRs
    model@genes <- append(model@genes, list(character(0L)))
    model@gprRules <- append(model@gprRules, list(""))

    # constraints
    model@constraints@coeff <- cbind(model@constraints@coeff,
                                   rep(0, constraint_num(model)))

    # subsystems
    model@subSys <- rbind(model@subSys,
                        rep(0, ncol(model@subSys)))
  }

  # Add/update metabolites if necessary
  model <- addMetabolite(model = model, id = met, name = metName, comp = metComp,
                         chemicalFormula = metChemicalFormula, charge = metCharge)
  indMs <- match(met, model@met_id)

  # Add/update subsystems  if necessary
  if(length(subSystem) > 0) {
    # model <- addSubSystem(model, ...) TODO
    indSubSys <- match(subSystem, model@subSys_id)
  }

  #----------------------------#
  # Add/update reaction values #
  #----------------------------#

  # S and obj.-coeff.
  if(!is.na(Scoef[1]))
    model@S[matrix(c(rep(indR,length(indMs)), indMs),ncol = 2)] <- Scoef
  if(!is.na(obj))
    model@obj_coef[indR] <- obj

  # reaction slots
  if(!is.na(reactName))
    model@react_name[indR] <- reactName
  if(!is.na(annotation))
    model@react_attr[indR]$annotation <- annotation

  # bounds
  model@lowbnd[indR] <- lb
  model@uppbnd[indR] <- ub

  # subsys
  if(length(subSystem) > 0) {
    model@S[matrix(c(indSubSys, rep(indR,length(indSubSys))),ncol = 2)] <- TRUE
  }

  return(model)
}

#' Add metabolites or update their data
#'
#' The functions allows you to add one or more metabolites to a model. When
#' providing the ID of an already existing metabolite, you can use this function
#' to update metabolite information.
#'
#' @param model Model of class \link{modelorg}
#'
#' @export
addMetabolite <- function(model, id, name = NA, comp = NA, chemicalFormula = NA,
                          charge = NA, annotation = NA) {
  # TODO
}
