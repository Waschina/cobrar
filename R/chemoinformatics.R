#' Count elements in formulas
#'
#' Counts the number of elements in molecules based on their chemical formulas.
#'
#' @param formulas Character vector of chemical formulas.
#' @param multiplier Multiplier for the element coutns in each molecule
#'
#' @returns A numeric matrix with named columns for individual elements. Rows
#' contain the counts of elements for each formula in the same order as formulas
#' were provided to the function.
#'
#' @details
#' Formulas may not have element numbers containing decimal points.
#'
#'
#' @examples
#' countElements(c("C6H12O6","C48H72CoN11O8","HCOOH"))
#'
#' @importFrom stats aggregate
#'
#' @export
countElements <- function(formulas, multiplier = 1) {

  if(length(multiplier) == 1) {
    multiplier <- rep(multiplier, length(formulas))
  }

  if(length(multiplier) != length(formulas))
    stop("Number of provided formulas does not match mumber of multipliers.")

  formulas <- ifelse(formulas == "",NA_character_,formulas)

  rgx <- "[A-Z][a-z]*[0-9]*"
  elects <- gregexec(rgx, formulas)
  elects <- lapply(1:length(formulas), function(i) {
    tmp <- list(ects = elects[[i]],
                formula = formulas[i])
    return(tmp)
  })

  out <- lapply(elects, function(x) {
    epos <- x$ects
    elen <- attributes(x$ects)$match.length
    eend <- epos + elen - 1

    tmpcts <- c()
    for(j in 1:length(epos)) {
      res <- substr(x$formula,epos[j],eend[j])
      tmpcts <- c(tmpcts, res)
    }

    if(length(tmpcts) == 0)
      return(numeric(0))

    ects <- as.numeric(gsub("^[A-z]*","",tmpcts))
    ects <- ifelse(is.na(ects),1,ects)
    eche <- gsub("[0-9]*$","",tmpcts)

    acts <- aggregate(ects, by = list(eche), sum)

    res <- acts[,"x"]
    names(res) <- acts[,1]
    if(length(res) == 0)
      res <- numeric(0)
    return(res)
  })

  uniqElements <- unique(unlist(lapply(out, names)))

  eleMarkup <- matrix(0, nrow = length(formulas), ncol = length(uniqElements),
                      dimnames = list(NULL,uniqElements))

  for(i in 1:length(out)) {
    eleMarkup[i,names(out[[i]])] <- out[[i]]
  }

  eleMarkup <- eleMarkup * multiplier

  return(eleMarkup)
}

tmp <- countElements(c("C6H12O6","FeO",NA,"CH"), multiplier = c(1,-1,3,2))


