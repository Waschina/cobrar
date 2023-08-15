setClass("LPproblem",
         slots = c(
           ptr = "externalptr",
           solver = "character",
           method = "character"
         )
)


# Set generic functions, where each solver gets its own version

setGeneric("loadLPprob", function(lp, ...) {
  standardGeneric("loadLPprob")
})

setGeneric("setObjDirection", function(lp, ...) {
  standardGeneric("setObjDirection")
})

setGeneric("addCols", function(lp, ...) {
  standardGeneric("addCols")
})

setGeneric("addRows", function(lp, ...) {
  standardGeneric("addRows")
})

setGeneric("loadMatrix", function(lp, ...) {
  standardGeneric("loadMatrix")
})

setGeneric("setColsBndsObjCoefs", function(lp, ...) {
  standardGeneric("setColsBndsObjCoefs")
})

setGeneric("setColsKind", function(lp, ...) {
  standardGeneric("setColsKind")
})

setGeneric("setRowsBnds", function(lp, ...) {
  standardGeneric("setRowsBnds")
})

setGeneric("solveLp", function(lp, ...) {
  standardGeneric("solveLp")
})

setGeneric("getObjValue", function(lp, ...) {
  standardGeneric("getObjValue")
})

setGeneric("getSolStat", function(lp, ...) {
  standardGeneric("getSolStat")
})

setGeneric("getColsPrimal", function(lp, ...) {
  standardGeneric("getColsPrimal")
})

setGeneric("getRedCosts", function(lp, ...) {
  standardGeneric("getRedCosts")
})

setGeneric("returnCode", function(lp, ...) {
  standardGeneric("returnCode")
})

setGeneric("statusCode", function(lp, ...) {
  standardGeneric("statusCode")
})

setGeneric("addSingleConstraint", function(lp, ...) {
  standardGeneric("addSingleConstraint")
})

# setGeneric("addConstraints", function(lp, ...) {
#   standardGeneric("addConstraints")
# })

setGeneric("fvaJob", function(lp, ...) {
  standardGeneric("fvaJob")
})
