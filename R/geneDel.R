#' Identify reactions affected by gene knockouts
#'
#' This function identifies reactions, which cannot be catalyzed anymore when
#' a specified set of genes is deleted from a model.
#'
#' @param model Model of class \link{ModelOrg}
#' @param gene Character or numeric vector providing the IDs or indices of genes
#' to be deleted from 'model'.
#' @param single If FALSE (default), the effect of simultaneous gene deletions
#' of all genes in `gene` are assumed. If TRUE, results for single gene
#' deletions are returned.
#'
#' @returns Character vector with reactions IDs if `single` is FALSE, an a list
#' of character vectors if `single` is TRUE. In the latter case, the i-th
#' element in the returned list corresponds to the gene deletion results of the
#' i-th gene of the input parameter.
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#'
#' # Identify reactions that would be lost after multiple gene deletions
#' geneDel(mod, gene = c("b3916","b1723","b4025"))
#'
#' # Identify reactions that would be lost after single gene deletions
#' geneDel(mod, gene = c("b3916","b1723","b4025"), single = TRUE)
#'
#' @export
geneDel <- function(model, gene, single = FALSE) {

  # check if provides gene values are valid
  if(any(!checkGeneId(model, gene = gene))) {
    stop("check argument gene.")
  }

  if(is.numeric(gene))
    gene <- model@allGenes[gene]

  KOgenes <- gene

  # Identify reactions associated to at least one of the genes of interest
  if(!single) {
    roi_pos <- which(unlist(lapply(model@genes, FUN = function(x) any(gene %in% x))))
    roi_pos <- list(roi_pos)
    KOgenes <- list(KOgenes)
  } else {
    roi_pos <- lapply(gene, FUN = function(x) {
      which(unlist(lapply(model@genes, FUN = function(y) any(y == x))))
    })
    KOgenes <- as.list(KOgenes)
  }

  res <- lapply(1:length(roi_pos), FUN = function(i) {
    roi_pos_i <- roi_pos[[i]]
    if(length(roi_pos_i) == 0)
      return(character(0L))
    gprL <- lapply(roi_pos_i,
                   function(i) c(gpr = model@gprRules[i], genes = model@genes[i]))

    # can each reaction be catalysed after gene KO(s)?
    reactCata <- unlist(lapply(gprL, function(gpr_i) {
      if(is.null(gpr_i$gpr) || gpr_i$gpr == "")
        return(TRUE)
      x <- gpr_i$genes
      x <- !(x %in% KOgenes[[i]] | is.na(x))
      return(eval(parse(text=gpr_i$gpr)))
    }))

    return(model@react_id[roi_pos_i[!reactCata]])
  })

  if(!single) {
    res <- unlist(res)
  }

  return(res)
}
