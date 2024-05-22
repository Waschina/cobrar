BMterms <- c(
  "Growth", # carveme models
  "bio1", # gapseq/modelseed
  "BIOMASS.*" # bigg(?)
)

#'
#' Guess biomass reaction(s)
#'
#' @param model Model of class \link{ModelOrg}
#'
#' @docType methods
#' @rdname guessBMReaction-methods
#' @export
setGeneric("guessBMReaction", valueClass = "character", function(model) {
  standardGeneric("guessBMReaction")
})

#' @rdname guessBMReaction-methods
#' @aliases guessBMReaction,ModelOrg
setMethod("guessBMReaction", signature(model = "ModelOrg"),
          function(model) {
            termregex <- paste0("^",BMterms,"$", collapse = "|")
            ind <- grep(termregex, model@react_id)
            return(model@react_id[ind])
          }
)

#' @rdname guessBMReaction-methods
#' @aliases guessBMReaction,ModelComm
setMethod("guessBMReaction", signature(model = "ModelComm"),
          function(model) {
            termregex <- paste0("^M[0-9]+_",BMterms,"$", collapse = "|")
            ind <- grep(termregex, model@react_id)
            return(model@react_id[ind])
          }
)
