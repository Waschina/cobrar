#' Reads a sybil model file and constructs an object of cobrar's class
#' 'ModelOrg'
#'
#' @param file_path Path to sybil model(s) saved as RDS file.
#'
#' @returns A \link{ModelOrg-class} object or a list with \link{ModelOrg-class}
#' objects.
#'
#' @import Matrix
#' @family Model import and export helpers
#' @export
readSybilmod <- function(file_path) {
  sybildoc.lst  <- readRDS(normalizePath(file_path))

  if( !is.list(sybildoc.lst) ) sybildoc.lst <- list(sybildoc.lst)

  models.lst <- vector(mode="list", length=length(sybildoc.lst))
  names(models.lst) <- names(sybildoc.lst)
  for(i in seq_along(sybildoc.lst)){
    sybildoc <- sybildoc.lst[[i]]
    constraints <- new("Constraints",
                       coeff = as(Matrix(nrow = 0, ncol = ncol(sybildoc@S), sparse = TRUE),
                                  "dMatrix"),
                       lb = numeric(0),
                       ub = numeric(0),
                       rtype = character(0))

    mod <- new("ModelOrg",
               mod_id = sybildoc@mod_id,
               mod_desc = sybildoc@mod_desc,
               mod_name = sybildoc@mod_name,
               mod_compart = sybildoc@mod_compart,
               mod_compart_name = sybildoc@mod_compart,
               mod_attr = sybildoc@mod_attr,
               mod_notes = paste0("<notes>\n  <html xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>",
                                  sybildoc@mod_desc,
                                  "</p>\n  </html>\n</notes>"),
               S = sybildoc@S,
               obj_coef = sybildoc@obj_coef,
               subSys = sybildoc@subSys,
               subSys_id = colnames(sybildoc@subSys),
               subSys_name = colnames(sybildoc@subSys),
               constraints = constraints,

               met_id = sybildoc@met_id,
               met_name = sybildoc@met_name,
               met_comp = sybildoc@mod_compart[sybildoc@met_comp],
               met_attr = sybildoc@met_attr,

               react_id = sybildoc@react_id,
               react_name = sybildoc@react_name,
               react_comp = "",
               lowbnd = sybildoc@lowbnd,
               uppbnd = sybildoc@uppbnd,
               react_attr = sybildoc@react_attr,

               gprRules = sybildoc@gprRules,
               genes = sybildoc@genes,
               allGenes = sybildoc@allGenes,
               genes_attr = data.frame(name = sybildoc@allGenes,
                                       CVTerms = "", SBOTerm = ""))

    # add new columns to data table
    mod@mod_attr$CVTerms <- ""
    mod@mod_attr$SBOTerm <- ""
    mod@met_attr$CVTerms <- ""
    mod@met_attr$SBOTerm <- ""
    mod@react_attr$CVTerms <- ""
    mod@react_attr$SBOTerm <- ""

    models.lst[[i]] <- mod
  }

  if( length(models.lst) == 1 ){
    return(models.lst[[1]])
  }else{
    return(models.lst)
  }
}
