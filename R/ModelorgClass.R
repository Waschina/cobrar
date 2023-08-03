#' Stucture of modelorg Class
#'
#' This class represents a model organization with various attributes related
#' to central model structures, metabolites, reactions, and genes.
#'
#' @slot mod_id A character vector representing the model identifier.
#' @slot mod_desc A character vector describing the model.
#' @slot mod_name A character vector containing the model name.
#' @slot mod_compart A character vector indicating the model compartment.
#' @slot mod_compartName A character vector with the name of the model compartment.
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
#' @export
setClass("modelorg",

         slots = c(
           # central model structures
           mod_id = "character",
           mod_desc = "character",
           mod_name = "character",
           mod_compart = "character",
           mod_compartName = "character",
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
