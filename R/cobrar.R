## usethis namespace: start
#' @useDynLib cobrar, .registration = TRUE
#' @importFrom Rcpp evalCpp
## usethis namespace: end
NULL

.onAttach <- function(libname, pkgname) {

  packageStartupMessage("cobrar. Let's do this.")
}
