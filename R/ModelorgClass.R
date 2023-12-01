#' Structure of modelorg Class
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
#' @slot constraints An object of class \link{Constraints} which specifies
#' constraints in a model in addition to stationarity and individual flux
#' bounds.
#'
#' @slot met_id A character vector representing metabolite identifiers.
#' @slot met_name A character vector with metabolite names.
#' @slot met_comp A character vector indicating metabolite compartments.
#' @slot met_attr A data.frame that enables the storage of additional data for
#' metabolites. Only specific columns are exported to SBML files. See
#' \link{writeSBMLmod} for details.
#'
#' @slot react_id A character vector representing reaction identifiers.
#' @slot react_name A character vector with reaction names.
#' @slot react_comp A character vector indicating reaction compartments.
#' @slot lowbnd A character vector containing lower bounds for reactions.
#' @slot uppbnd A character vector containing upper bounds for reactions.
#' @slot react_attr A data.frame that enables the storage of additional data for
#' reactions. Only specific columns are exported to SBML files. See
#' \link{writeSBMLmod} for details.
#'
#' @slot gprRules A character vector with Gene-Protein-Reaction association rules
#' (with gene product indices corresponding to the order in slot 'genes').
#' @slot genes A list of character vectors. Each vector contains the IDs of
#' gene products associated to the respective reaction.
#' @slot allGenes A character vector with all gene identifiers.
#' @slot genes_attr A data.frame that enables the storage of additional data
#' (e.g., name and CVTerms) for genes/gene products. Only specific columns are
#' exported to SBML files. See \link{writeSBMLmod} for details.
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
           constraints = "Constraints",

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
           genes_attr = "data.frame"
         ))

#------------------------------------------------------------------------------#
# setters and getters                                                          #
#------------------------------------------------------------------------------#

#' Number of reactions
#'
#' Get the total number of reactions of a model
#'
#' @param model Model of class \link{modelorg}
#'
#' @export
setGeneric("react_num", valueClass = "numeric", function(model) {
  standardGeneric("react_num")
})
setMethod("react_num", signature(model = "modelorg"),
          function(model) {
            return(length(model@react_id))
          }
)

#' Index of reaction(s)
#'
#' Returns the index(es) of specific reaction(s).
#'
#' @param model Model of class \link{modelorg}
#' @param react Character vector with reaction IDs or Integer vector providing
#' indexes.
#'
#' @details
#' Returns NA for reaction IDs not part of the model or the index is larger than
#' the number of reactions in the model.
#'
#' @export
setGeneric("react_pos", valueClass = "numeric", function(model, react) {
  standardGeneric("react_pos")
})
setMethod("react_pos", signature(model = "modelorg", react = "character"),
          function(model, react) {
            return(match(react, model@react_id))
          }
)
setMethod("react_pos", signature(model = "modelorg", react = "numeric"),
          function(model, react) {
            return(ifelse(react <= react_num(model),react, NA_integer_))
          }
)
setMethod("react_pos", signature(model = "modelorg", react = "missing"),
          function(model, react) {
            return(NA_integer_)
          }
)


#' Number of metabolites
#'
#' Get the total number of metabolites of a model
#'
#' @param model Model of class \link{modelorg}
#'
#' @export
setGeneric("met_num", valueClass = "numeric", function(model) {
  standardGeneric("met_num")
})
setMethod("met_num", signature(model = "modelorg"),
          function(model) {
            return(length(model@met_id))
          }
)

#' Index of metabolite(s)
#'
#' Returns the index(es) of specific metabolite(s).
#'
#' @param model Model of class \link{modelorg}
#' @param met Character vector with metabolite IDs or Integer vector providing
#' indexes.
#'
#' @details
#' Returns NA for metabolite IDs not part of the model or if the index is larger
#' than the number of metabolites in the model.
#'
#' @export
setGeneric("met_pos", valueClass = "numeric", function(model, met) {
  standardGeneric("met_pos")
})
setMethod("met_pos", signature(model = "modelorg", met = "character"),
          function(model, met) {
            return(match(met, model@met_id))
          }
)
setMethod("met_pos", signature(model = "modelorg", met = "numeric"),
          function(model, met) {
            return(ifelse((met<=met_num(model)), met, NA_integer_))
          }
)
setMethod("met_pos", signature(model = "modelorg", met = "missing"),
          function(model, met) {
            return(NA_integer_)
          }
)

#' Number of genes
#'
#' Get the total number of genes of a model
#'
#' @param model Model of class \link{modelorg}
#'
#' @export
setGeneric("gene_num", valueClass = "numeric", function(model) {
  standardGeneric("gene_num")
})
setMethod("gene_num", signature(model = "modelorg"),
          function(model) {
            return(length(model@allGenes))
          }
)


#' Index of gene(s)
#'
#' Returns the index(es) of specific gene(s).
#'
#' @param model Model of class \link{modelorg}
#' @param gene Character vector with gene IDs or Integer vector providing
#' indexes.
#'
#' @details
#' Returns NA for gene IDs not part of the model or if the index is larger
#' than the number of genes in the model.
#'
#' @export
setGeneric("gene_pos", valueClass = "numeric", function(model, gene) {
  standardGeneric("gene_pos")
})
setMethod("gene_pos", signature(model = "modelorg", gene = "character"),
          function(model, gene) {
            return(match(gene, model@allGenes))
          }
)
setMethod("gene_pos", signature(model = "modelorg", gene = "numeric"),
          function(model, gene) {
            return(ifelse(gene<=gene_num(model),gene,NA_integer_))
          }
)
setMethod("gene_pos", signature(model = "modelorg", gene = "missing"),
          function(model, gene) {
            return(NA_integer_)
          }
)

#' Number of compartments
#'
#' Get the total number of compartments of a model
#'
#' @param model Model of class \link{modelorg}
#'
#' @export
setGeneric("comp_num", valueClass = "numeric", function(model) {
  standardGeneric("comp_num")
})
setMethod("comp_num", signature(model = "modelorg"),
          function(model) {
            return(length(model@mod_compart))
          }
)

#' Index of compartment(s)
#'
#' Returns the index(es) of specific compartment(s).
#'
#' @param model Model of class \link{modelorg}
#' @param comp Character vector with compartment IDs or Integer vector providing
#' indexes.
#'
#' @details
#' Returns NA for compartment IDs not part of the model or if the index is larger
#' than the number of compartments in the model.
#'
#' @export
setGeneric("comp_pos", valueClass = "numeric", function(model, comp) {
  standardGeneric("comp_pos")
})
setMethod("comp_pos", signature(model = "modelorg", comp = "character"),
          function(model, comp) {
            return(match(comp, model@mod_compart))
          }
)
setMethod("comp_pos", signature(model = "modelorg", comp = "numeric"),
          function(model, comp) {
            return(ifelse(comp<=comp_num(model),comp,NA_integer_))
          }
)
setMethod("comp_pos", signature(model = "modelorg", comp = "missing"),
          function(model, comp) {
            return(NA_integer_)
          }
)
setMethod("comp_pos", signature(model = "modelorg", comp = "logical"),
          function(model, comp) {
            return(rep(NA_integer_, length(comp)))
          }
)

#' Number of constraints
#'
#' Get the total number of constraints of a model
#'
#' @param model Model of class \link{modelorg}
#'
#' @export
setGeneric("constraint_num", valueClass = "numeric", function(model) {
  standardGeneric("constraint_num")
})
setMethod("constraint_num", signature(model = "modelorg"),
          function(model) {
            return(nrow(model@constraints@coeff))
          }
)

#' Number of subsystems
#'
#' Get the total number of subsystems of a model
#'
#' @param model Model of class \link{modelorg}
#'
#' @export
setGeneric("subsys_num", valueClass = "numeric", function(model) {
  standardGeneric("subsys_num")
})
setMethod("subsys_num", signature(model = "modelorg"),
          function(model) {
            return(length(model@subSys_id))
          }
)

#' Index of subsystem(s)
#'
#' Returns the index(es) of specific subsystem(s).
#'
#' @param model Model of class \link{modelorg}
#' @param subsys Character vector with subsystem IDs or Integer vector providing
#' indexes.
#'
#' @details
#' Returns NA for subsystem IDs not part of the model or if the index is larger
#' than the number of subsystems in the model.
#'
#' @export
setGeneric("subsys_pos", valueClass = "numeric", function(model, subsys) {
  standardGeneric("subsys_pos")
})
setMethod("subsys_pos", signature(model = "modelorg", subsys = "character"),
          function(model, subsys) {
            return(match(subsys, model@subSys_id))
          }
)
setMethod("subsys_pos", signature(model = "modelorg", subsys = "numeric"),
          function(model, subsys) {
            return(ifelse(subsys<=subsys_num(model),subsys,NA_integer_))
          }
)
setMethod("subsys_pos", signature(model = "modelorg", subsys = "missing"),
          function(model, subsys) {
            return(NA_integer_)
          }
)
setMethod("subsys_pos", signature(model = "modelorg", subsys = "logical"),
          function(model, subsys) {
            return(rep(NA_integer_, length(subsys)))
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
            cat("number of user constraints:", constraint_num(object), "\n")
            if(constraint_num(object) > 0 && constraint_num(object) <= 10) {
              for(i in 1:constraint_num(object)) {
                cat("                           ",constraint2string(object,i),"\n")
              }
            }
            cat("number of subsystems:      ", subsys_num(object), "\n")
            cat("\n")
            cat("objective function:        ", printObjFunc(object), "\n")
          }
)

setGeneric("constraint2string" ,valueClass = "character", function(object, ind, ...) {
  standardGeneric("constraint2string")
})
setMethod("constraint2string", signature(object = "modelorg", ind = "numeric"),
          function(object, ind, digits = 5) {
            nz <- which(object@constraints@coeff[ind,] != 0)
            cnz <- c()
            for(i in 1:length(nz)) {
              cnz[i] <- paste0(ifelse(sign(object@constraints@coeff[ind,nz[i]])==1,"+","-"),
                               round(abs(object@constraints@coeff[ind,nz[i]]), digits = digits)," ",
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
            ccstr <- sapply(1:constraint_num(object), function(i) constraint2string(object , i, digits = Inf))

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
