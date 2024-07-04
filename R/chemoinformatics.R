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

#' Data frame of elements and their average weights
#'
#' Gets a data.frame with all elements, their symbol, name, atomic number and
#' their abridged standard atomic weight according to the International Union of
#' Pure and Applied Chemistry (IUPAC).
#'
#' @references
#' Prohaska T, Irrgeher J, Benefield J, Böhlke JK, Chesson LA, Coplen TB, et al.
#' Standard atomic weights of the elements 2021 (IUPAC Technical Report). Vol.
#' 94, Pure and Applied Chemistry. 2022. p. 573–600.
#' http://dx.doi.org/10.1515/pac-2019-0603
#'
#' @export
elements <- function() {
  return(.COBRARenv$elements)
}


#' Calculate molar mass of molecules
#'
#' Calculates the average molar mass of compounds based on their chemical
#' formulas.
#'
#' @param formulas Character vector of chemical formulas.
#'
#' @examples
#' mass(c("C6H12O6","C48H72CoN11O8","HCOOH"))
#'
#' @export
mass <- function(formulas) {
  makeup <- countElements(formulas)

  ele <- elements()
  eleMasses <- ele$AvgWeight
  names(eleMasses) <- ele$Symbol

  # are there any unknown/unstable elements in the makeup matrix?
  indNoElement <- which(!(colnames(makeup) %in% ele$Symbol))
  indUnstableElement <- which(colnames(makeup) %in% ele$Symbol[is.na(ele$AvgWeight)])

  if(length(indNoElement)>0)
    warning(paste0("Some element symbols in the formulas are not known: ",
                   paste(colnames(makeup)[indNoElement],collapse = ", "),
                   ". Omitting those symbols in mass calculation."))

  if(length(indUnstableElement)>0)
    warning(paste0("Some elements in the formulas do not have stable nuclides: ",
                   paste(colnames(makeup)[indUnstableElement],collapse = ", "),
                   ". Omitting those elements in mass calculation."))

  colkeep <- !((1:ncol(makeup)) %in% c(indNoElement,indUnstableElement))
  makeup <- makeup[, colkeep, drop = FALSE]
  masses <- (makeup %*% eleMasses[colnames(makeup)])[,1]

  return(masses)
}

