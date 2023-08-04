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
#'
#' @slot met_id A character vector representing metabolite identifiers.
#' @slot met_name A character vector with metabolite names.
#' @slot met_comp A character vector indicating metabolite compartments.
#' @slot met_attr A character vector with additional metabolite attributes.
#'
#' @slot react_id A character vector representing reaction identifiers.
#' @slot react_name A character vector with reaction names.
#' @slot lowbnd A character vector containing lower bounds for reactions.
#' @slot uppbnd A character vector containing upper bounds for reactions.
#'
#' @slot gprRules A character vector with Gene-Protein-Reaction association rules.
#' @slot genes A list of gene-related information.
#' @slot allGenes A character vector with all gene identifiers.
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

           # metabolites,
           met_id = "character",
           met_name = "character",
           met_comp = "integer",
           met_attr = "data.frame",

           # reactions
           react_id = "character",
           react_name = "character",
           react_comp = "integer",
           lowbnd = "numeric",
           uppbnd = "numeric",
           react_attr = "data.frame",

           # genes
           gprRules = "character",
           genes = "list",
           allGenes = "character"
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
            cat("model ID:              ", object@mod_id, "\n")
            cat("model name:            ", object@mod_name, "\n")
            cat("number of compartments:", length(object@mod_compart), "\n")
            for(i in 1:length(object@mod_compart)) {
              cat("                       ", object@mod_compart[i], " (",
                  object@mod_compart_name[i], ")\n")
            }
            cat("number of reactions:   ", react_num(object), "\n")
            cat("number of metabolites: ", met_num(object), "\n")
            if (length(object@allGenes) > 0) {
              cat("number of unique genes:", length(object@allGenes), "\n")
            }
            cat("objective function:    ", printObjFunc(object), "\n")
          }
)


