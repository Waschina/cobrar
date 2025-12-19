#' Find exchange reactions
#'
#' Finds all exchange reactions within a (community) model and report them in a
#' data frame.
#'
#' @param model Model of class \link{ModelOrg} or \link{ModelComm}
#'
#' @returns If 'model' is of class \link{ModelOrg}, a data.frame is returned
#' stating all exchange reaction IDs, their index in the reaction list, the
#' respective metabolite id, name, and index in the metabolite list. If the
#' model is of class \link{ModelComm}, the resulting data.frame contains the
#' community exchange reactions and the organism-specific exchange reaction. The
#' latter are exchange reactions, that connect extracellular metabolites of the
#' organism metabolic network models with the shared extracellular space.
#'
#' @details
#' In cobrar, an exchange reaction is an unbalanced reaction involving a single
#' metabolite that appears exclusively on one side of the reaction equation
#' (either as a substrate or as a product, but not both). This structure
#' represents the import or export of that metabolite between the system and its
#' environment.
#'
#'
#' @docType methods
#'
#' @rdname findExchReact-methods
#'
#' @export
setGeneric("findExchReact" ,valueClass = "data.frame", function(model) {
  standardGeneric("findExchReact")
})
#' @rdname findExchReact-methods
#' @aliases findExchReact,ModelOrg
setMethod("findExchReact", signature(model = "ModelOrg"),
          function(model) {
            ex_pos <- which(diff(model@S@p) == 1)
            ex_id <- model@react_id[ex_pos]

            cpd_pos <- apply(model@S[,ex_pos],2,function(x) which(x!=0))
            cpd_id <- model@met_id[cpd_pos]
            cpd_name <- model@met_name[cpd_pos]

            return(data.frame(react_id = ex_id,
                              react_pos = ex_pos,
                              met_id = cpd_id,
                              met_name = cpd_name,
                              met_pos = cpd_pos,
                              lb = model@lowbnd[ex_pos],
                              ub = model@uppbnd[ex_pos]))
          })
#' @rdname findExchReact-methods
#' @aliases findExchReact,ModelComm
setMethod("findExchReact", signature(model = "ModelComm"),
          function(model) {
            ex_pos <- which(diff(model@S@p) == 1)
            ex_id <- model@react_id[ex_pos]

            # mex_id <- model@react_id[grep("^M[0-9]+_EX_",model@react_id)]
            # mex_pos <- react_pos(model, mex_id)
            emets <- which(model@met_comp == "e")
            Stmp <- model@S[-emets,]
            mex_pos <- which(diff(Stmp@p) == 1)
            mex_id <- model@react_id[mex_pos]

            cpd_posEX <- apply(model@S[,ex_pos],2,function(x) which(x!=0))
            cpd_posMEX <- apply(model@S[,mex_pos],2,function(x) which(x!=0 & model@met_comp != "e"))
            cpd_pos <- c(cpd_posEX, cpd_posMEX)
            cpd_id <- model@met_id[cpd_pos]
            cpd_name <- model@met_name[cpd_pos]

            return(data.frame(exchange_type = c(rep("Community",length(ex_id)),
                                                rep("Organism",length(mex_id))),
                              organism_id = c(rep(NA_character_,length(ex_id)),
                                              sub("(^M[0-9]+)_.*$","\\1",mex_id)),
                              react_id = c(ex_id, mex_id),
                              react_pos = c(ex_pos, mex_pos),
                              met_id = cpd_id,
                              met_name = cpd_name,
                              met_pos = cpd_pos,
                              lb = model@lowbnd[c(ex_pos,mex_pos)],
                              ub = model@uppbnd[c(ex_pos,mex_pos)]))
          })
