#' Stucture of UserConstraints Class
#'
#' This class represents user constraints that can be added to a model of class
#' \link{modelorg} in addition to the stationarity constraint (\eqn{S v = 0})
#' and flux bounds.
#'
#' @slot coeff A sparse numeric matrix of \link[Matrix]{dgCMatrix-class}
#' representing the coefficients for each reaction in the model. Each row
#' denotes a user constraint, each column a reaction in the model in the same
#' order as in slor "S" in the corresponding \link{modelorg} object.
#' @slot lb Numeric vector provididing the lower bound for each constraint.
#' @slot ub Numeric vector provididing the lower bound for each constraint.
#' @slot rtype Character vector stating the constraint type. See details.
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
#' @aliases UserConstraints
#'
#' @exportClass UserConstraints
setClass("UserConstraints",
         slots = c(
           coeff = "dgCMatrix",
           lb = "numeric",
           ub = "numeric",
           rtype = "character"
         )
)

#' Stucture of modelorg Class
#'
#' This class represents a model organization with various attributes related
#' to central model structures, metabolites, reactions, and genes.
#'
#' @slot mod_id A character vector representing the model identifier.
#' @slot mod_desc A character vector describing the model.
#' @slot mod_name A character vector containing the model name.
#' @slot mod_compart A character vector indicating the model compartment.
#' @slot mod_compart_name A character vector with the name of the model compartment.
#' @slot mod_attr A data frame with additional model attributes.
#' @slot S A sparse numeric matrix of \link[Matrix]{dgCMatrix-class} representing the Stoichiometric matrix.
#' @slot obj_coef A numeric vector containing coefficients for the objective function.
#' @slot subSys A sparse Boolean matrix of \link[Matrix]{lgCMatrix-class} defining subsystems.
#' @slot subSys_id A character vector representing subsystem identifiers.
#' @slot subSys_name A character vector containing the subsystem names.
#' @slot constraints An object of class \link{UserConstraints} which specifies
#' constraints in a model in addition to stationarity and individual flux
#' bounds.
#'
#' @slot met_id A character vector representing metabolite identifiers.
#' @slot met_name A character vector with metabolite names.
#' @slot met_comp A character vector indicating metabolite compartments.
#' @slot met_attr A character vector with additional metabolite attributes.
#'
#' @slot react_id A character vector representing reaction identifiers.
#' @slot react_name A character vector with reaction names.
#' @slot react_comp A character vector indicating reaction compartments.
#' @slot lowbnd A character vector containing lower bounds for reactions.
#' @slot uppbnd A character vector containing upper bounds for reactions.
#'
#' @slot gprRules A character vector with Gene-Protein-Reaction association rules
#' (with gene product indexes corresponding to the order in slot 'genes').
#' @slot genes A list of character vectors. Each vector contains the IDs of
#' gene products associated to the respective reaction.
#' @slot allGenes A character vector with all gene identifiers.
#' @slot allGenes_name A character vector with all gene names.
#'
#' @aliases modelorg
#'
#' @exportClass modelorg
setClass("modelorg",

         slots = c(
           # central model structures
           mod_id = "character",
           mod_desc = "character",
           mod_name = "character",
           mod_compart = "character",
           mod_compart_name = "character",
           mod_attr = "data.frame",
           mod_notes = "character",
           S = "dgCMatrix",
           obj_coef = "numeric",
           subSys = "lgCMatrix",
           subSys_id = "character",
           subSys_name = "character",
           constraints = "UserConstraints",

           # metabolites,
           met_id = "character",
           met_name = "character",
           met_comp = "character",
           met_attr = "data.frame",

           # reactions
           react_id = "character",
           react_name = "character",
           react_comp = "character",
           lowbnd = "numeric",
           uppbnd = "numeric",
           react_attr = "data.frame",

           # genes
           gprRules = "character",
           genes = "list",
           allGenes = "character",
           allGenes_name = "character"
         ))

#------------------------------------------------------------------------------#
# setters and getters                                                          #
#------------------------------------------------------------------------------#

# number of reactions
setGeneric("react_num", valueClass = "numeric", function(object) {
  standardGeneric("react_num")
})
setMethod("react_num", signature(object = "modelorg"),
          function(object) {
            return(length(object@react_id))
          }
)

# index of reaction(s)
setGeneric("react_pos", valueClass = "numeric", function(object, react) {
  standardGeneric("react_pos")
})
setMethod("react_pos", signature(object = "modelorg", react = "character"),
          function(object, react) {
            return(match(react, object@react_id))
          }
)
setMethod("react_pos", signature(object = "modelorg", react = "numeric"),
          function(object, react) {
            return(react)
          }
)

# number of metabolites
setGeneric("met_num", valueClass = "numeric", function(object) {
  standardGeneric("met_num")
})
setMethod("met_num", signature(object = "modelorg"),
          function(object) {
            return(length(object@met_id))
          }
)

# index of metabolite(s)
setGeneric("met_pos", valueClass = "numeric", function(object, met) {
  standardGeneric("met_pos")
})
setMethod("met_pos", signature(object = "modelorg", met = "character"),
          function(object, met) {
            return(match(met, object@met_id))
          }
)
setMethod("met_pos", signature(object = "modelorg", met = "numeric"),
          function(object, met) {
            return(met)
          }
)


# number of genes
setGeneric("gene_num", valueClass = "numeric", function(object) {
  standardGeneric("gene_num")
})
setMethod("gene_num", signature(object = "modelorg"),
          function(object) {
            return(length(object@allGenes))
          }
)


# index of gene(s)
setGeneric("gene_pos", valueClass = "numeric", function(object, gene) {
  standardGeneric("gene_pos")
})
setMethod("gene_pos", signature(object = "modelorg", gene = "character"),
          function(object, gene) {
            return(match(gene, object@allGenes))
          }
)
setMethod("gene_pos", signature(object = "modelorg", gene = "numeric"),
          function(object, gene) {
            return(gene)
          }
)


# number of user constraints
setGeneric("constraints_num", valueClass = "numeric", function(object) {
  standardGeneric("constraints_num")
})
setMethod("constraints_num", signature(object = "modelorg"),
          function(object) {
            return(nrow(object@constraints@coeff))
          }
)

#------------------------------------------------------------------------------#
# Miscellaneous                                                                #
#------------------------------------------------------------------------------#

setGeneric("printObjFunc", valueClass = "character", function(object) {
  standardGeneric("printObjFunc")
})
setMethod("printObjFunc", signature(object = "modelorg"),
          function(object) {
            cInd <- object@obj_coef != 0

            # check if there is an objective function
            if (sum(cInd) == 0) {
              of <- "no objective function"
            }
            else {
              obj <- gsub("^([^-])", "+\\1",
                          object@obj_coef[cInd], perl = TRUE)
              of  <- paste(paste(obj, object@react_id[cInd]),
                           collapse = " ")
            }

            return(of)
          }
)

#' Print a short summary of a metabolic network model
#'
#' Displays a few key properties of a metabolic network model of class
#' \link{modelorg}.
#'
#' @param object S4-object of class \link{modelorg}.
#'
#' @export
setMethod("show", signature(object = "modelorg"),
          function(object) {
            cat("model ID:                  ", object@mod_id, "\n")
            cat("model name:                ", object@mod_name, "\n")
            cat("number of compartments:    ", length(object@mod_compart), "\n")
            for(i in 1:length(object@mod_compart)) {
              cat("                           ", object@mod_compart[i], " (",
                  object@mod_compart_name[i], ")\n")
            }
            cat("number of reactions:       ", react_num(object), "\n")
            cat("number of metabolites:     ", met_num(object), "\n")
            cat("number of unique genes:    ", gene_num(object), "\n")
            cat("number of user constraints:", constraints_num(object), "\n")
            if(constraints_num(object) > 0 && constraints_num(object) <= 10) {
              for(i in 1:constraints_num(object)) {
                cat("                           ",constraint2string(object,i),"\n")
              }
            }
            cat("\n")
            cat("objective function:        ", printObjFunc(object), "\n")
          }
)

setGeneric("constraint2string" ,valueClass = "character", function(object, ind) {
  standardGeneric("constraint2string")
})
setMethod("constraint2string", signature(object = "modelorg", ind = "integer"),
          function(object, ind) {
            nz <- which(object@constraints@coeff[ind,] != 0)
            cnz <- c()
            for(i in 1:length(nz)) {
              cnz[i] <- paste0(ifelse(sign(object@constraints@coeff[ind,nz[i]])==1,"+","-"),
                               round(abs(object@constraints@coeff[ind,nz[i]]), digits = 5)," ",
                               object@react_id[nz[i]])
            }
            mid <- paste(cnz, collapse = " ")
            lhs <- switch(object@constraints@rtype[ind],
              "F" = "-Inf < ",
              "L" = paste0(object@constraints@lb[ind]," <= "),
              "U" = "-Inf < ",
              "D" = paste0(object@constraints@lb[ind]," <= "),
              "E" = paste0(object@constraints@lb[ind]," == ")
            )

            rhs <- switch(object@constraints@rtype[ind],
                          "F" = " < Inf",
                          "U" = paste0(" <= ", object@constraints@ub[ind]),
                          "L" = " < Inf",
                          "D" = paste0(" <= ", object@constraints@ub[ind]),
                          "E" = paste0(" == ", object@constraints@ub[ind])
            )

            return(paste0(lhs, mid, rhs))
          }
)


# remove duplicate constraints
setGeneric("rmDuplicateConstraints", valueClass = "modelorg", function(object) {
  standardGeneric("rmDuplicateConstraints")
})
setMethod("rmDuplicateConstraints", signature(object = "modelorg"),
          function(object) {
            ccstr <- sapply(1:constraints_num(object), function(i) constraint2string(object , i))

            indrm <- which(duplicated(ccstr))

            if(length(indrm)>0) {
              warning("Duplicate user constraints. Retaining only unique constraints.")
              object@constraints@coeff <- object@constraints@coeff[-indrm,, drop = FALSE]
              object@constraints@lb <- object@constraints@lb[-indrm]
              object@constraints@ub <- object@constraints@ub[-indrm]
              object@constraints@rtype <- object@constraints@rtype[-indrm]
            }
            return(object)
          }
)
