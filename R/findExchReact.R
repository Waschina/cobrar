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
#' @export
setGeneric("findExchReact" ,valueClass = "data.frame", function(model) {
  standardGeneric("findExchReact")
})
setMethod("findExchReact", signature(model = "ModelOrg"),
          function(model) {
            ex_id <- model@react_id[grep("^EX_",model@react_id)]
            ex_pos <- react_pos(model, ex_id)

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
setMethod("findExchReact", signature(model = "ModelComm"),
          function(model) {
            ex_id <- model@react_id[grep("^EX_",model@react_id)]
            ex_pos <- react_pos(model, ex_id)

            mex_id <- model@react_id[grep("^M[0-9]+_EX_",model@react_id)]
            mex_pos <- react_pos(model, mex_id)

            cpd_posEX <- apply(model@S[,ex_pos],2,function(x) which(x!=0))
            cpd_posMEX <- apply(model@S[grep("^M[0-9]+",model@met_id),mex_pos],2,function(x) which(x!=0))
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
