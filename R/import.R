#' Reads an sybil model file and constructs an object of class 'modelorg'
#'
#' @param file_path Path to sybil model(s) saved as RDS file.
#'
#' @returns A \link{modelorg-class} object or a list with \link{modelorg-class} objects.
#'
#' @import Matrix
#' @export
readSybilmod <- function(file_path) {
  sybildoc.lst  <- readRDS(normalizePath(file_path))

  if( !is.list(sybildoc.lst) ) sybildoc.lst <- list(sybildoc.lst)

  models.lst <- list()
  for(sybildoc in sybildoc.lst){
      constraints <- new("Constraints",
                         coeff = as(Matrix(nrow = 0, ncol = ncol(sybildoc@S), sparse = TRUE),
                                    "dMatrix"),
                         lb = numeric(0),
                         ub = numeric(0),
                         rtype = character(0))

    mod <- new("modelorg",
                mod_id = sybildoc@mod_id,
                mod_desc = sybildoc@mod_desc,
                mod_name = sybildoc@mod_name,
                mod_compart = sybildoc@mod_compart,
                mod_compart_name = "",
                mod_attr = sybildoc@mod_attr,
                mod_notes = "",
                S = sybildoc@S,
                obj_coef = sybildoc@obj_coef,
                subSys = sybildoc@subSys,
                subSys_id = "",
                subSys_name = "",
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
                genes_attr = data.frame())
 
    models.lst <- append(models.lst, mod) 
  }
  
  if( length(models.lst) == 1 ){
    return(models.lst[[1]])
  }else{
    return(models.lst)
  }
}
