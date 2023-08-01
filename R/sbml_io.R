#' Reads an SBML file and constructs an object of class 'modelorg'
#'
#' @param filename Path to SBML file.
#'
#' @export
readSBMLmod <- function(filename) {
  sbmldoc <- readSBMLfile(filename)

  return(sbmldoc)
}
