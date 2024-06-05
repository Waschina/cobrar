#' @exportClass LPproblem
setClass("LPproblem",
         slots = c(
           ptr = "externalptr",
           solver = "character",
           method = "character"
         )
)


# Set generic functions, where each solver gets its own version
#' @export
setGeneric("loadLPprob", function(lp, ...) {
  standardGeneric("loadLPprob")
})
#' @export
setGeneric("setObjDirection", function(lp, ...) {
  standardGeneric("setObjDirection")
})
#' @export
setGeneric("addCols", function(lp, ...) {
  standardGeneric("addCols")
})
#' @export
setGeneric("addRows", function(lp, ...) {
  standardGeneric("addRows")
})
#' @export
setGeneric("loadMatrix", function(lp, ...) {
  standardGeneric("loadMatrix")
})
#' @export
setGeneric("setColsBndsObjCoefs", function(lp, ...) {
  standardGeneric("setColsBndsObjCoefs")
})
#' @export
setGeneric("setColsKind", function(lp, ...) {
  standardGeneric("setColsKind")
})
#' @export
setGeneric("setRowsBnds", function(lp, ...) {
  standardGeneric("setRowsBnds")
})
#' @export
setGeneric("solveLp", function(lp, ...) {
  standardGeneric("solveLp")
})
#' @export
setGeneric("getObjValue", function(lp, ...) {
  standardGeneric("getObjValue")
})
#' @export
setGeneric("getSolStat", function(lp, ...) {
  standardGeneric("getSolStat")
})
#' @export
setGeneric("getColsPrimal", function(lp, ...) {
  standardGeneric("getColsPrimal")
})
#' @export
setGeneric("getRedCosts", function(lp, ...) {
  standardGeneric("getRedCosts")
})
#' @export
setGeneric("returnCode", function(lp, ...) {
  standardGeneric("returnCode")
})
#' @export
setGeneric("statusCode", function(lp, ...) {
  standardGeneric("statusCode")
})
#' @export
setGeneric("addSingleConstraint", function(lp, ...) {
  standardGeneric("addSingleConstraint")
})

# setGeneric("addConstraints", function(lp, ...) {
#   standardGeneric("addConstraints")
# })
#' @export
setGeneric("fvaJob", function(lp, ...) {
  standardGeneric("fvaJob")
})
