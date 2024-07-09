#' Reads an SBML file and constructs an object of class 'ModelOrg'
#'
#' @param file_path Path to SBML file.
#'
#' @returns A \link{ModelOrg-class} object.
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#' mod
#'
#' @import Matrix
#' @export
readSBMLmod <- function(file_path) {
  is_compressed <- FALSE

  file_path <- normalizePath(file_path)

  if(grepl("\\.gz$",file_path)) {
    is_compressed <- TRUE
    fileconn <- gzfile(file_path)
    tmpfile <- tempfile(fileext = ".xml")
    writeLines(readLines(fileconn), tmpfile)
    close(fileconn)
    sbmldoc  <- readSBMLfile(tmpfile)
  } else {
    sbmldoc  <- readSBMLfile(file_path)
  }

  # Pointers to SBML document and model (libSBML objects)
  modelPtr <- getModelObj(sbmldoc)

  #---------------#
  # Model content #
  #---------------#

  # Model fields
  mod_id <- getModelId(modelPtr)
  mod_name <- getModelName(modelPtr)
  if(is.na(mod_name))
    mod_name <- mod_id
  S <- getStoichiometricMatrix(modelPtr)
  mod_compartments <- getModelCompartments(modelPtr)
  mod_compartments$name <- ifelse(is.na(mod_compartments$name),
                                  mod_compartments$id,
                                  mod_compartments$name)
  mod_cvterms <- paste(getModelCVTerms(modelPtr), collapse = ";")
  mod_notes <- getModelNotes(modelPtr)
  obj <- getObjectiveFunction(modelPtr)
  constraints <- new("Constraints",
                     coeff = as(Matrix(nrow = 0, ncol = ncol(S), sparse = TRUE),
                                "dMatrix"),
                     lb = numeric(0),
                     ub = numeric(0),
                     rtype = character(0))
  mod_sbo <- getModelSBOTerm(modelPtr)

  # Subsystems
  subSys <- getSubsystems(modelPtr); colnames(subSys$subSys) <- subSys$subSys_ids

  # Reactions
  react_id <- getReactionIds(modelPtr)
  react_name <- getReactionNames(modelPtr)
  react_bnds <- getReactionFluxBounds(modelPtr)
  react_comp <- getReactionCompartment(modelPtr)
  react_cvterms <- getReactionCVTerms(modelPtr)
  react_cvterms <- lapply(react_cvterms, function(x) {
    paste(x, collapse = ";")
  })
  react_sboterms <- getReactionSBOTerms(modelPtr)

  # Metabolites
  met_id <- getMetaboliteIds(modelPtr)
  met_name <- getMetaboliteNames(modelPtr)
  met_attr <- getMetaboliteAnnotation(modelPtr)
  met_cvterms <- getMetaboliteCVTerms(modelPtr)
  met_cvterms <- lapply(met_cvterms, function(x) {
    paste(x, collapse = ";")
  })
  met_attr$CVTerms <- unlist(met_cvterms)
  met_comp <- getMetaboliteCompartments(modelPtr)
  met_attr$SBOTerm <- getMetaboliteSBOTerms(modelPtr)


  # Genes
  allGeneProducts <- getGeneProducts(modelPtr)
  gpr <- getGPRs(modelPtr)
  gpr_cvterms <- getGeneProductCVTerms(modelPtr)
  gpr_cvterms <- lapply(gpr_cvterms, function(x) {
    paste(x, collapse = ";")
  })
  gpr_sbo <- getGeneProductSBOTerms(modelPtr)

  if(is_compressed)
    file.remove(tmpfile)

  return(
    new("ModelOrg",
        mod_id = mod_id,
        mod_desc = mod_id,
        mod_name = mod_name,
        mod_compart = mod_compartments$id,
        mod_compart_name = mod_compartments$name,
        mod_attr = data.frame(CVTerms = mod_cvterms, SBOTerm = mod_sbo),
        mod_notes = mod_notes,
        S = S,
        obj_coef = obj$coeff,
        obj_dir = ifelse(!(obj$dir %in% c("minimize","maximize")),
                         "maximize", obj$dir),
        subSys = as(subSys$subSys, "lMatrix"),
        subSys_id = subSys$subSys_ids,
        subSys_name = subSys$subSys_names,
        constraints = constraints,

        met_id = gsub("^M_","",met_id),
        met_name = met_name,
        met_comp = ifelse(met_comp == "", NA_character_, met_comp),
        met_attr = met_attr,

        react_id = gsub("^R_","",react_id),
        react_name = react_name,
        react_comp = ifelse(react_comp == "", NA_character_, react_comp),
        lowbnd = react_bnds$lower_bound,
        uppbnd = react_bnds$upper_bound,
        react_attr = data.frame(CVTerms = unlist(react_cvterms),
                                SBOTerm = react_sboterms),

        gprRules = gpr$rules,
        genes = lapply(gpr$genes, function(x) gsub("^G_","",x)),
        allGenes = gsub("^G_","",allGeneProducts$ID),
        genes_attr = data.frame(name = allGeneProducts$name,
                                CVTerms = unlist(gpr_cvterms),
                                SBOTerm = gpr_sbo)
    )
  )
}


# Small helpter function to transform an SBO Term to it's integer als ID
sboterm2int <- function(sbo) {
  return(as.numeric(gsub("SBO:|SBO_","",sbo)))
}

#' Exports a Metabolic Network in SBML Format
#'
#' Export a constraint-based metabolic network model from a S4 object of class
#' \link{ModelOrg} to a SBML file.
#'
#' @param model Model of class \link{ModelOrg}
#' @param file_path SBML file name for exporting the model. Default is the
#' model's ID with ".xml" suffix.
#'
#' @details
#' Exported SBML files are of level 3, version 2. FBC-package version 2.
#'
#' What content from the data.frames `react_attr`, `met_attr`, and `mod_attr` is
#' exported to SBML files? Currently only the columns named "CVTerms".
#'
#'
#' @returns TRUE if file export was successful.
#'
#' @export
writeSBMLmod <- function(model, file_path = NULL) {

  compress <- FALSE
  if(!is.null(file_path) && grepl("\\.gz$",file_path)) {
    compress <- TRUE
    out_file <- file_path
    file_path <- gsub("\\.gz$","", file_path)
  }

  if(is.null(file_path))
    file_path <- paste0(model@mod_id, ".xml")

  file_path <- path.expand(file_path)

  # small corrections before export
  if(!all(grepl("^R_",model@react_id)))
    model@react_id <- paste0("R_",model@react_id)
  if(!all(grepl("^M_",model@met_id)))
    model@met_id <- paste0("M_",model@met_id)
  if(!all(grepl("^G_",model@allGenes))) {
    model@allGenes <- paste0("G_",model@allGenes)
    model@genes <- lapply(model@genes, function(x) paste0("G_",x))
  }
  colnames(model@subSys) <- NULL
  model@subSys_id <- gsub("-","_",model@subSys_id)
  if(!all(grepl("^subsys_",model@subSys_id)))
    model@subSys_id <- paste0("subsys_",model@subSys_id)
  model@met_id <- sub("\\[(.*)\\]$", "_\\1", model@met_id)

  # libSBML seems not to allow "." in gene IDs... replacing them here
  if(any(grepl("\\.",model@allGenes))) {
    model@allGenes <- gsub("\\.","_",model@allGenes)
    model@genes <- lapply(model@genes, FUN = function(x) gsub("\\.","_",x))
    # warning("Some gene IDs contain dots ('.'). Replacing them with underscores ('_').")
  }


  if(is.na(model@mod_attr$CVTerms[1]))
    model@mod_attr$CVTerms[1] <- ""

  if(react_num(model) > 0)
    model@react_attr$CVTerms <- ifelse(is.na(model@react_attr$CVTerms),
                                       "",model@react_attr$CVTerms)
  if(met_num(model) > 0)
    model@met_attr$CVTerms <- ifelse(is.na(model@met_attr$CVTerms),
                                     "",model@met_attr$CVTerms)
  if(gene_num(model) > 0)
    model@genes_attr$CVTerms <- ifelse(is.na(model@genes_attr$CVTerms),
                                       "",model@genes_attr$CVTerms)

  # Stoichiometry lists
  lReaMets <- apply(model@S, 2, FUN = function(x) model@met_id[which(abs(x)>0)])
  lReaStoich <- apply(model@S, 2, FUN = function(x) x[which(abs(x)>0)])

  # bound groups
  bndgrp <- data.frame(id = model@react_id, lb = model@lowbnd, ub = model@uppbnd,
                       lb.term = rep(NA_character_, react_num(model)),
                       ub.term = rep(NA_character_, react_num(model)))
  if(nrow(bndgrp)>0){
    bndgrp$lb.term <- paste0(bndgrp$id,"_lb")
    bndgrp$ub.term <- paste0(bndgrp$id,"_ub")
    bndgrp$lb.term <- ifelse(bndgrp$lb == 0,"default_0", bndgrp$lb.term)
    bndgrp$ub.term <- ifelse(bndgrp$ub == 0,"default_0", bndgrp$ub.term)
    bndgrp$lb.term <- ifelse(bndgrp$lb == -COBRAR_SETTINGS("MAXIMUM"),"default_lb", bndgrp$lb.term)
    bndgrp$ub.term <- ifelse(bndgrp$ub == COBRAR_SETTINGS("MAXIMUM"),"default_ub", bndgrp$ub.term)
  }
  bndgrp_para <- data.frame(bnd = c(bndgrp$lb.term, bndgrp$ub.term),
                            val = c(bndgrp$lb, bndgrp$ub))
  bndgrp_para <- bndgrp_para[!duplicated(bndgrp_para$bnd),]
  bndgrp_para$SBO <- ifelse(grepl("^default_",bndgrp_para$bnd),626,625)
  bndgrp_para <- bndgrp_para[order(bndgrp_para$bnd),]

  # gpr string for libSBML
  gpr <- character(0L)
  if(react_num(model)>0) {
    for(i in 1:length(model@react_id)) {
      x <- model@gprRules[[i]]
      y <- model@genes[[i]]
      y <- ifelse(y == "G_NA","false",y)

      if(length(y) == 0) {
        gpr <- append(gpr, "")
      } else {
        for(j in 1:length(y)) {
          x <- gsub(paste0("x[",j,"]"),
                    y[j], x, fixed = TRUE)
        }
        gpr <- append(gpr, x)
      }
    }
    gpr <- gsub("\\&", "and", gpr)
    gpr <- gsub("\\|", "or", gpr)
  }
  #print(gpr)

  # Let's export
  out <- writeSBML(
    file_path = file_path,

    # Model fields
    mod_id = model@mod_id,
    mod_name = model@mod_name,
    mod_desc = model@mod_desc,
    mod_cvterms = unlist(strsplit(model@mod_attr$CVTerms[1], ";")),
    mod_notes = model@mod_notes,
    mod_sbo = sboterm2int(model@mod_attr$SBOTerm[1]),

    # Compartments
    comp_id = model@mod_compart,
    comp_name = model@mod_compart_name,

    # Species
    met_id = model@met_id,
    met_name = model@met_name,
    met_charge = model@met_attr$charge,
    met_formula = model@met_attr$chemicalFormula,
    met_comp = model@met_comp,
    met_cvterms = strsplit(model@met_attr$CVTerms, ";"),
    met_sbo = sboterm2int(model@met_attr$SBOTerm),

    # Subsystems
    subsys = apply(model@subSys, 2, function(x) which(x)-1),
    subsys_id = model@subSys_id,
    subsys_name = model@subSys_name,

    # Genes
    gene_id = model@allGenes,
    gene_name = model@genes_attr$name,
    gene_cvterms = strsplit(model@genes_attr$CVTerms, ";"),
    gene_sbo = sboterm2int(model@genes_attr$SBOTerm),

    # Parameters (bound groups)
    param_id = bndgrp_para$bnd,
    param_val = bndgrp_para$val,
    param_sbo = bndgrp_para$SBO,

    # Reactions and Stoichiometry
    react_id = model@react_id,
    react_name = model@react_name,
    Scoeff = lReaStoich,
    react_mets = lReaMets,
    react_lb = bndgrp$lb.term,
    react_ub = bndgrp$ub.term,
    react_rev = ifelse(bndgrp$lb < 0 & bndgrp$ub > 0, TRUE, FALSE),
    react_cvterms = strsplit(model@react_attr$CVTerms, ";"),
    react_sbo = sboterm2int(model@react_attr$SBOTerm),
    gpr = gpr,

    # Objective
    obj_coef = model@obj_coef,
    obj_dir = model@obj_dir
  )

  if(compress) {
    input_conn <- file(file_path, "r")
    output_conn <- gzfile(out_file, "w")
    file_content <- readLines(input_conn)
    writeLines(file_content, output_conn)
    close(input_conn)
    close(output_conn)
    file.remove(file_path)
  }

  return(out)
}


#' Reads a sybil model file and constructs an object of cobrar's class
#' 'ModelOrg'
#'
#' @param file_path Path to sybil model(s) saved as RDS file.
#'
#' @returns A \link{ModelOrg-class} object or a list with \link{ModelOrg-class}
#' objects.
#'
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
               obj_dir = "maximize",
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
