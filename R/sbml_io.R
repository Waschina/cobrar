#' Reads an SBML file and constructs an object of class 'modelorg'
#'
#' @param file_path Path to SBML file.
#'
#' @returns A \link{modelorg-class} object.
#'
#' @examples
#' library(cobrar)
#'
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#'
#' @import Matrix
#' @export
readSBMLmod <- function(file_path) {

  # Pointers to SBML document and model (libSBML objects)
  sbmldoc  <- readSBMLfile(file_path)
  modelPtr <- getModelObj(sbmldoc)

  #---------------#
  # Model content #
  #---------------#

  # Model fields
  mod_id <- getModelId(modelPtr)
  mod_name <- getModelName(modelPtr)
  S <- getStoichiometricMatrix(modelPtr)
  mod_compartments <- getModelCompartments(modelPtr)
  mod_annotation <- getModelAnnotation(modelPtr)
  mod_notes <- getModelNotes(modelPtr)
  obj_coeff <- getObjectiveFunction(modelPtr)
  subSys <- getSubsystems(modelPtr); colnames(subSys$subSys) <- subSys$subSys_ids

  # Reactions
  react_id <- getReactionIds(modelPtr)
  react_name <- getReactionNames(modelPtr)
  react_bnds <- getReactionFluxBounds(modelPtr)
  react_anno <- getReactionAnnotation(modelPtr)

  # Metabolites
  met_id <- getMetaboliteIds(modelPtr)
  met_name <- getMetaboliteNames(modelPtr)
  met_attr <- getMetaboliteAnnotation(modelPtr)

  # Genes
  # TODO

  # Clean up (frees memory)
  rm(sbmldoc, modelPtr)

  return(
    new("modelorg",
        mod_id = mod_id,
        mod_desc = mod_id,
        mod_name = mod_name,
        mod_compart = mod_compartments$id,
        mod_compartName = mod_compartments$name,
        mod_attr = data.frame(annotation = mod_annotation),
        mod_notes = mod_notes,
        S = S,
        obj_coef = obj_coeff,
        subSys = as(subSys$subSys, "lMatrix"),
        subSys_id = subSys$subSys_ids,
        subSys_name = subSys$subSys_names,

        met_id = met_id,
        met_name = met_name,
        met_comp = NA_integer_, # TODO
        met_attr = met_attr,

        react_id = react_id,
        react_name = react_name,
        react_comp = NA_integer_, # TODO
        lowbnd = react_bnds$lower_bound,
        uppbnd = react_bnds$upper_bound,
        react_attr = data.frame(annotation = react_anno),

        gprRules = NA_character_, # TODO
        genes = list(), # TODO
        allGenes = NA_character_ # TODO
    )
  )
}
