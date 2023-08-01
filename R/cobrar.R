## usethis namespace: start
#' @useDynLib cobrar, .registration = TRUE
## usethis namespace: end
NULL

#' @importFrom utils installed.packages
.onAttach <- function(libname, pkgname) {

  packageStartupMessage("cobrar. Let's do this.")
}
