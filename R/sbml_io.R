#' Reads an SBML file and constructs an object of class 'modelorg'
#'
#' @param filename Path to SBML file.
#'
#' @export
readSBMLmod <- function(file_path) {

  # Pointers to SBML document and model (libSBML objects)
  sbmldoc  <- readSBMLfile(file_path)
  modelPtr <- getModelObj(sbmldoc)

  #---------------#
  # Model content #
  #---------------#

  # Stoichiometric matrix
  S <- getStoichiometricMatrix(modelPtr)

  # Reactions
  rxn_id <- getReactionIds(modelPtr)
  rxn_names <- getReactionNames(modelPtr)
  rxn_bnds <- getReactionFluxBounds(modelPtr)

  return(list(
    pointer_sbml = sbmldoc,
    pointer_model = modelPtr,
    S = S,
    rxn_id = rxn_id,
    rxn_names = rxn_names,
    rxn_bnds = rxn_bnds
  ))
}
