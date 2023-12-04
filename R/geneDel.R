#' Identify reactions affected by gene knockouts
#'
#' This function identifies reactions, which cannot be catalyzed anymore when
#' a specified set of genes is deleted from a model.
#'
#' @param model Model of class \link{modelorg}
#' @param gene Character or numeric vector providing the IDs or indices of genes
#' to be deleted from 'model'.
#'
#' @return Character vector with reactions IDs.
#'
#' @export
geneDel <- function(model, gene) {

  # check if provides gene values are valid
  if(any(!checkGeneId(model, gene = gene))) {
    stop("check argument gene.")
  }

  KOgenes <- gene

  if(is.numeric(gene))
    KOgenes <- model@allGenes[gene]

  gprL <- lapply(1:react_num(model),
                 function(i) c(gpr = model@gprRules[i], genes = model@genes[i]))

  # can each reaction be catalysed after gene KO(s)?
  reactCata <- unlist(lapply(gprL, function(gpr_i) {
    if(is.null(gpr_i$gpr) || gpr_i$gpr == "")
      return(TRUE)
    x <- gpr_i$genes
    # return(x)
    x <- !(x %in% KOgenes)
    return(eval(parse(text=gpr_i$gpr)))
  }))

  return(model@react_id[!reactCata])
}
