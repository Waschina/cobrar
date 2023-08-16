## usethis namespace: start
#' @useDynLib cobrar, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom methods .valueClassTest as is new
## usethis namespace: end
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("cobrar uses...\n",
                        " - libSBML (v. ", getSBMLVersion(),")\n",
                        " - glpk (v. ",getGLPKVersion(),")")
}

.COBRARenv <- new.env()

.onLoad <- function(lib, pkg) {

  # central settings
  .COBRARenv$settings <- list(
    MODELORG_VERSION = "3.0",
    SOLVER           = "glpk",
    METHOD           = "simplex",
    TOLERANCE        = 1E-6,
    MAXIMUM          = 1000,
    ALGORITHM        = "fba",
    OPT_DIRECTION    = "max",
    USE_NAMES        = FALSE,
    PATH_TO_MODEL    = ".",
    SOLVER_CTRL_PARM = as.data.frame(NA)
  )

  # solvers
  .COBRARenv$solvers <- c("glpk") # hopefully IBM's cplex soon


  # methods
  .COBRARenv$solverMethods <- list(
    glpk = c("simplex", "interior", "exact", "mip")
  )

  # default parameters
  .COBRARenv$solverCtrlParm <- list(

    glpk = list(simplex  = as.data.frame(NA),
                interior = as.data.frame(NA),
                exact    = as.data.frame(NA),
                mip      = as.data.frame(NA)
    )
  )

  # algorithms simulation genetic perturbations
  .COBRARenv$algorithm[["pert"]] <- c("fba")

}

#' Set and get central cobrar parameters
#'
#' Manage a set of default parameter settings for cobrar.
#'
#' @param parm A character string giving the name of the parameter to set.
#' @param value Choose a value to over-write the current parameter.
#' @param ... TBD.
#'
#' @export
COBRAR_SETTINGS <- function(parm, value, ...) {

  if ( (missing(parm)) && (missing(value)) ) {
    return(.COBRARenv$settings)
  }

  if (missing(value)) {
    if (!parm %in% names(.COBRARenv$settings)) {
      stop("unknown parameter ", sQuote(parm))
    }
    return(.COBRARenv$settings[[parm]])
  }

  if ( (length(parm) != 1) ||
       ( (length(value) != 1) && (! (parm == "SOLVER_CTRL_PARM") ) ) ) {
    stop("arguments 'parm' and 'value' must have a length of 1.")
  }

  switch(parm,
         "MODELORG_VERSION" = {
           stop("this value must not be set by the user!")
         },

         "SOLVER" = {
           if(value != "glpk")
             warning("Currently, only glpk is supported as solver.")
         },

         "METHOD" = {
           if(value %in% .COBRARenv$solverMethods[[.COBRARenv$settings[["SOLVER"]]]]) {
             .COBRARenv$settings[["METHOD"]] <- value
           } else {
             warning(paste0(
               "Method ", value, " not available for solver ", .COBRARenv$settings[["SOLVER"]],". ",
               "Staying with current method: ", .COBRARenv$settings[["METHOD"]]
             ))
           }

           # chmet <- checkDefaultMethod(solver = COBRAR_SETTINGS("SOLVER"),
           #                             method = value,
           #                             probType = "", ...)
           # .COBRARenv$settings[["SOLVER"]]           <- chmet$sol
           # .COBRARenv$settings[["METHOD"]]           <- chmet$met
           # .COBRARenv$settings[["SOLVER_CTRL_PARM"]] <- chmet$parm
         },

         "TOLERANCE" = {
           .COBRARenv$settings[["TOLERANCE"]] <- as.numeric(value)
         },

         "MAXIMUM" = {
           .COBRARenv$settings[["MAXIMUM"]] <- as.numeric(value)
         },

         "ALGORITHM" = {
           .COBRARenv$settings[["ALGORITHM"]] <- as.character(value)
         },

         "OPT_DIRECTION" = {
           if ( (value == "max") || (value == "min") ) {
             .COBRARenv$settings[["OPT_DIRECTION"]] <- as.character(value)
           }
           else {
             stop("OPT_DIRECTION can be either 'max' or 'min'")
           }
         },

         "USE_NAMES" = {
           .COBRARenv$settings[["USE_NAMES"]] <- as.logical(value)
         },

         "PATH_TO_MODEL" = {
           if (file.exists(value)) {
             .COBRARenv$settings[["PATH_TO_MODEL"]] <- as.character(value)
           }
           else {
             stop("directory ", sQuote(value), " does not exist")
           }
         },

         "SOLVER_CTRL_PARM" = {
           if ( (is.data.frame(value)) || (is.list(value)) ) {
             if ("NA" %in% names(COBRAR_SETTINGS("SOLVER_CTRL_PARM"))) {
               .COBRARenv$settings[["SOLVER_CTRL_PARM"]] <- value
             }
             else {
               pn <- names(value)
               for (i in seq(along = value)) {
                 .COBRARenv$settings[["SOLVER_CTRL_PARM"]][[pn[i]]] <- value[[pn[i]]]
               }
             }
           }
           else {
             stop("SOLVER_CTRL_PARM must be data.frame or list")
           }
         },

         {
           stop("unknown parameter: ", sQuote(parm))
         }
  )
}
