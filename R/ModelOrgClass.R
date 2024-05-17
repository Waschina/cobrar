#' Structure of ModelOrg Class
#'
#' This class represents a model organism with various attributes related
#' to central model structures, metabolites, reactions, and genes.
#'
#' @slot mod_id A character vector representing the model identifier.
#' @slot mod_desc A character vector describing the model.
#' @slot mod_name A character vector containing the model name.
#' @slot mod_compart A character vector indicating the model compartment.
#' @slot mod_compart_name A character vector with the name of the model compartment.
#' @slot mod_attr A data frame with additional model attributes.
#' @slot mod_notes A character string that can contain an XML block with
#' additional information about the model.
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
#' @aliases ModelOrg
#'
#' @family Object classes
#' @exportClass ModelOrg
setClass("ModelOrg",

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
         ),
         prototype = list(
           # model slots
           mod_id = "cobrar-model",
           mod_desc = NA_character_,
           mod_name = "cobrar-model",
           mod_compart = character(0L),
           mod_compart_name = character(0L),
           mod_attr = data.frame(CVTerms = "", SBOTerm = "SBO:0000624"),
           mod_notes = NA_character_,
           S = as(as(as(Matrix(nrow = 0, ncol = 0, sparse = TRUE),
                        "dMatrix"),
                     "generalMatrix"),
                  "CsparseMatrix"),
           obj_coef = numeric(0L),
           subSys = as(as(as(Matrix(nrow = 0, ncol = 0, sparse = TRUE),
                             "lMatrix"),
                          "generalMatrix"),
                       "CsparseMatrix"),
           subSys_id = character(0L),
           subSys_name = character(0L),
           constraints = new("Constraints",
                             coeff = as(as(as(Matrix(nrow = 0, ncol = 0, sparse = TRUE),
                                        "dMatrix"),"generalMatrix"),"CsparseMatrix"),
                             lb = numeric(0),
                             ub = numeric(0),
                             rtype = character(0)),

           # Metabolites
           met_id = character(0L),
           met_name = character(0L),
           met_comp = character(0L),
           met_attr = data.frame(chemicalFormula = character(0L),
                                 charge = numeric(0L),
                                 CVTerms = character(0L),
                                 SBOTerm = character(0L)),

           # reactions
           react_id = character(0L),
           react_name = character(0L),
           react_comp = character(0L),
           lowbnd = numeric(0),
           uppbnd = numeric(0),
           react_attr = data.frame(CVTerms = character(0L),
                                   SBOTerm = character(0L)),

           # genes
           gprRules = character(0L),
           genes = list(),
           allGenes = character(0L),
           genes_attr = data.frame(name = character(0L),
                                   CVTerms = character(0L),
                                   SBOTerm = character(0L))
         ))


#------------------------------------------------------------------------------#
# setters and getters                                                          #
#------------------------------------------------------------------------------#

#' Number of reactions
#'
#' Get the total number of reactions of a model
#'
#' @param model Model of class \link{ModelOrg}
#'
#' @docType methods
#' @rdname react_num-methods
#' @family Model characteristics
#' @export
setGeneric("react_num", valueClass = "numeric", function(model) {
  standardGeneric("react_num")
})
#' @rdname react_num-methods
#' @aliases react_num,ModelOrg
setMethod("react_num", signature(model = "ModelOrg"),
          function(model) {
            return(length(model@react_id))
          }
)

#' Index of reaction(s)
#'
#' Returns the index(es) of specific reaction(s).
#'
#' @param model Model of class \link{ModelOrg}
#' @param react Character vector with reaction IDs or Integer vector providing
#' indexes.
#'
#' @details
#' Returns NA for reaction IDs not part of the model or the index is larger than
#' the number of reactions in the model.
#'
#' @docType methods
#' @rdname react_pos-methods
#' @family Model characteristics
#' @export
setGeneric("react_pos", valueClass = "numeric", function(model, react) {
  standardGeneric("react_pos")
})
#' @rdname react_pos-methods
#' @aliases react_pos,ModelOrg,character
setMethod("react_pos", signature(model = "ModelOrg", react = "character"),
          function(model, react) {
            return(match(react, model@react_id))
          }
)
#' @rdname react_pos-methods
#' @aliases react_pos,ModelOrg,numeric
setMethod("react_pos", signature(model = "ModelOrg", react = "numeric"),
          function(model, react) {
            return(ifelse(react <= react_num(model),react, NA_integer_))
          }
)
#' @rdname react_pos-methods
#' @aliases react_pos,ModelOrg,missing
setMethod("react_pos", signature(model = "ModelOrg", react = "missing"),
          function(model, react) {
            return(NA_integer_)
          }
)


#' Number of metabolites
#'
#' Get the total number of metabolites of a model
#'
#' @param model Model of class \link{ModelOrg}
#'
#' @docType methods
#' @rdname met_num-methods
#' @family Model characteristics
#' @export
setGeneric("met_num", valueClass = "numeric", function(model) {
  standardGeneric("met_num")
})
#' @rdname met_num-methods
#' @aliases met_num,ModelOrg
setMethod("met_num", signature(model = "ModelOrg"),
          function(model) {
            return(length(model@met_id))
          }
)

#' Index of metabolite(s)
#'
#' Returns the index(es) of specific metabolite(s).
#'
#' @param model Model of class \link{ModelOrg}
#' @param met Character vector with metabolite IDs or Integer vector providing
#' indexes.
#'
#' @details
#' Returns NA for metabolite IDs not part of the model or if the index is larger
#' than the number of metabolites in the model.
#'
#' @docType methods
#' @rdname met_pos-methods
#' @family Model characteristics
#' @export
setGeneric("met_pos", valueClass = "numeric", function(model, met) {
  standardGeneric("met_pos")
})
#' @rdname met_pos-methods
#' @aliases met_pos,ModelOrg,character
setMethod("met_pos", signature(model = "ModelOrg", met = "character"),
          function(model, met) {
            return(match(met, model@met_id))
          }
)
#' @rdname met_pos-methods
#' @aliases met_pos,ModelOrg,numeric
setMethod("met_pos", signature(model = "ModelOrg", met = "numeric"),
          function(model, met) {
            return(ifelse((met<=met_num(model)), met, NA_integer_))
          }
)
#' @rdname met_pos-methods
#' @aliases met_pos,ModelOrg,missing
setMethod("met_pos", signature(model = "ModelOrg", met = "missing"),
          function(model, met) {
            return(NA_integer_)
          }
)

#' Number of genes
#'
#' Get the total number of genes of a model
#'
#' @param model Model of class \link{ModelOrg}
#'
#' @docType methods
#' @rdname gene_num-methods
#' @family Model characteristics
#' @export
setGeneric("gene_num", valueClass = "numeric", function(model) {
  standardGeneric("gene_num")
})
#' @rdname gene_num-methods
#' @aliases gene_num,ModelOrg
setMethod("gene_num", signature(model = "ModelOrg"),
          function(model) {
            return(length(model@allGenes))
          }
)


#' Index of gene(s)
#'
#' Returns the index(es) of specific gene(s).
#'
#' @param model Model of class \link{ModelOrg}
#' @param gene Character vector with gene IDs or Integer vector providing
#' indexes.
#'
#' @details
#' Returns NA for gene IDs not part of the model or if the index is larger
#' than the number of genes in the model.
#'
#' @docType methods
#' @rdname gene_pos-methods
#' @family Model characteristics
#' @export
setGeneric("gene_pos", valueClass = "numeric", function(model, gene) {
  standardGeneric("gene_pos")
})
#' @rdname gene_pos-methods
#' @aliases gene_pos,ModelOrg,character
setMethod("gene_pos", signature(model = "ModelOrg", gene = "character"),
          function(model, gene) {
            return(match(gene, model@allGenes))
          }
)
#' @rdname gene_pos-methods
#' @aliases gene_pos,ModelOrg,numeric
setMethod("gene_pos", signature(model = "ModelOrg", gene = "numeric"),
          function(model, gene) {
            return(ifelse(gene<=gene_num(model),gene,NA_integer_))
          }
)
#' @rdname gene_pos-methods
#' @aliases gene_pos,ModelOrg,missing
setMethod("gene_pos", signature(model = "ModelOrg", gene = "missing"),
          function(model, gene) {
            return(NA_integer_)
          }
)

#' Number of compartments
#'
#' Get the total number of compartments of a model
#'
#' @param model Model of class \link{ModelOrg}
#'
#' @docType methods
#' @rdname comp_num-methods
#' @family Model characteristics
#' @export
setGeneric("comp_num", valueClass = "numeric", function(model) {
  standardGeneric("comp_num")
})
#' @rdname comp_num-methods
#' @aliases comp_num,ModelOrg
setMethod("comp_num", signature(model = "ModelOrg"),
          function(model) {
            return(length(model@mod_compart))
          }
)

#' Index of compartment(s)
#'
#' Returns the index(es) of specific compartment(s).
#'
#' @param model Model of class \link{ModelOrg}
#' @param comp Character vector with compartment IDs or Integer vector providing
#' indexes.
#'
#' @details
#' Returns NA for compartment IDs not part of the model or if the index is larger
#' than the number of compartments in the model.
#'
#' @docType methods
#' @rdname comp_pos-methods
#' @family Model characteristics
#' @export
setGeneric("comp_pos", valueClass = "numeric", function(model, comp) {
  standardGeneric("comp_pos")
})
#' @rdname comp_pos-methods
#' @aliases comp_pos,ModelOrg,character
setMethod("comp_pos", signature(model = "ModelOrg", comp = "character"),
          function(model, comp) {
            return(match(comp, model@mod_compart))
          }
)
#' @rdname comp_pos-methods
#' @aliases comp_pos,ModelOrg,numeric
setMethod("comp_pos", signature(model = "ModelOrg", comp = "numeric"),
          function(model, comp) {
            return(ifelse(comp<=comp_num(model),comp,NA_integer_))
          }
)
#' @rdname comp_pos-methods
#' @aliases comp_pos,ModelOrg,missing
setMethod("comp_pos", signature(model = "ModelOrg", comp = "missing"),
          function(model, comp) {
            return(NA_integer_)
          }
)
#' @rdname comp_pos-methods
#' @aliases comp_pos,ModelOrg,logical
setMethod("comp_pos", signature(model = "ModelOrg", comp = "logical"),
          function(model, comp) {
            return(rep(NA_integer_, length(comp)))
          }
)

#' Number of constraints
#'
#' Get the total number of constraints of a model
#'
#' @param model Model of class \link{ModelOrg}
#'
#' @docType methods
#' @rdname constraint_num-methods
#' @family Model characteristics
#' @export
setGeneric("constraint_num", valueClass = "numeric", function(model) {
  standardGeneric("constraint_num")
})
#' @rdname constraint_num-methods
#' @aliases constraint_num,ModelOrg
setMethod("constraint_num", signature(model = "ModelOrg"),
          function(model) {
            return(nrow(model@constraints@coeff))
          }
)

#' Number of subsystems
#'
#' Get the total number of subsystems of a model
#'
#' @param model Model of class \link{ModelOrg}
#'
#' @docType methods
#' @rdname subsys_num-methods
#' @family Model characteristics
#' @export
setGeneric("subsys_num", valueClass = "numeric", function(model) {
  standardGeneric("subsys_num")
})
#' @rdname subsys_num-methods
#' @aliases subsys_num,ModelOrg
setMethod("subsys_num", signature(model = "ModelOrg"),
          function(model) {
            return(length(model@subSys_id))
          }
)

#' Index of subsystem(s)
#'
#' Returns the index(es) of specific subsystem(s).
#'
#' @param model Model of class \link{ModelOrg}
#' @param subsys Character vector with subsystem IDs or Integer vector providing
#' indexes.
#'
#' @details
#' Returns NA for subsystem IDs not part of the model or if the index is larger
#' than the number of subsystems in the model.
#'
#' @docType methods
#' @rdname subsys_pos-methods
#' @family Model characteristics
#' @export
setGeneric("subsys_pos", valueClass = "numeric", function(model, subsys) {
  standardGeneric("subsys_pos")
})
#' @rdname subsys_pos-methods
#' @aliases subsys_pos,ModelOrg,character
setMethod("subsys_pos", signature(model = "ModelOrg", subsys = "character"),
          function(model, subsys) {
            return(match(subsys, model@subSys_id))
          }
)
#' @rdname subsys_pos-methods
#' @aliases subsys_pos,ModelOrg,numeric
setMethod("subsys_pos", signature(model = "ModelOrg", subsys = "numeric"),
          function(model, subsys) {
            return(ifelse(subsys<=subsys_num(model),subsys,NA_integer_))
          }
)
#' @rdname subsys_pos-methods
#' @aliases subsys_pos,ModelOrg,missing
setMethod("subsys_pos", signature(model = "ModelOrg", subsys = "missing"),
          function(model, subsys) {
            return(NA_integer_)
          }
)
#' @rdname subsys_pos-methods
#' @aliases subsys_pos,ModelOrg,logical
setMethod("subsys_pos", signature(model = "ModelOrg", subsys = "logical"),
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
setMethod("printObjFunc", signature(object = "ModelOrg"),
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
#' \link{ModelOrg}.
#'
#' @param object S4-object of class \link{ModelOrg}.
#'
#' @family Model characteristics
#' @export
setMethod("show", signature(object = "ModelOrg"),
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
setMethod("constraint2string", signature(object = "ModelOrg", ind = "numeric"),
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
setGeneric("rmDuplicateConstraints", valueClass = "ModelOrg", function(object) {
  standardGeneric("rmDuplicateConstraints")
})
setMethod("rmDuplicateConstraints", signature(object = "ModelOrg"),
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
