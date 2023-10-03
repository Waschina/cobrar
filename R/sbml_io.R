#' Reads an SBML file and constructs an object of class 'modelorg'
#'
#' @param file_path Path to SBML file.
#'
#' @returns A \link{modelorg-class} object.
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#' mod
#'
#' @import Matrix
#' @export
readSBMLmod <- function(file_path) {

  # Pointers to SBML document and model (libSBML objects)
  sbmldoc  <- readSBMLfile(normalizePath(file_path))
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
  mod_annotation <- getModelAnnotation(modelPtr)
  mod_notes <- getModelNotes(modelPtr)
  obj_coeff <- getObjectiveFunction(modelPtr)
  subSys <- getSubsystems(modelPtr); colnames(subSys$subSys) <- subSys$subSys_ids
  constraints <- new("Constraints",
                     coeff = as(Matrix(nrow = 0, ncol = ncol(S), sparse = TRUE),
                                "dMatrix"),
                     lb = numeric(0),
                     ub = numeric(0),
                     rtype = character(0))

  # Reactions
  react_id <- getReactionIds(modelPtr)
  react_name <- getReactionNames(modelPtr)
  react_bnds <- getReactionFluxBounds(modelPtr)
  # react_anno <- getReactionAnnotation(modelPtr)
  react_anno <- rep(NA_character_,ncol(S))
  react_comp <- getReactionCompartment(modelPtr)
  react_cvterms <- getReactionCVTerms(modelPtr)
  react_cvterms <- lapply(react_cvterms, function(x) {
    paste(x, collapse = ";")
  })

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


  # Genes
  allGeneProducts <- getGeneProducts(modelPtr)
  gpr <- getGPRs(modelPtr)

  return(
    new("modelorg",
        mod_id = mod_id,
        mod_desc = mod_id,
        mod_name = mod_name,
        mod_compart = mod_compartments$id,
        mod_compart_name = mod_compartments$name,
        mod_attr = data.frame(annotation = mod_annotation),
        mod_notes = mod_notes,
        S = S,
        obj_coef = obj_coeff,
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
        react_attr = data.frame(CVTerms = unlist(react_cvterms)),

        gprRules = gpr$rules,
        genes = lapply(gpr$genes, function(x) gsub("^G_","",x)),
        allGenes = gsub("^G_","",allGeneProducts$ID),
        allGenes_name = allGeneProducts$name
    )
  )
}

#' Exports a Metabolic Network in SBML Format
#'
#' Export a constraint-based metabolic network model from a S4 object of class
#' \link{modelorg} to a SBML file.
#'
#' @param model Model of class \link{modelorg}
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

  if(is.null(file_path))
    file_path <- paste0(model@mod_id, ".xml")

  # small corrections before export
  if(!all(grepl("^R_",model@react_id)))
    model@react_id <- paste0("R_",model@react_id)
  if(!all(grepl("^M_",model@met_id)))
    model@met_id <- paste0("M_",model@met_id)


  # Stoichiometry lists
  lReaMets <- apply(model@S, 2, FUN = function(x) model@met_id[which(abs(x)>0)])
  lReaStoich <- apply(model@S, 2, FUN = function(x) x[which(abs(x)>0)])

  # bound groups
  bndgrp <- data.frame(id = model@react_id, lb = model@lowbnd, ub = model@uppbnd)
  bndgrp$lb.term <- paste0(bndgrp$id,"_lb")
  bndgrp$ub.term <- paste0(bndgrp$id,"_ub")
  bndgrp$lb.term <- ifelse(bndgrp$lb == 0,"default_0", bndgrp$lb.term)
  bndgrp$ub.term <- ifelse(bndgrp$ub == 0,"default_0", bndgrp$ub.term)
  bndgrp$lb.term <- ifelse(bndgrp$lb == -COBRAR_SETTINGS("MAXIMUM"),"default_lb", bndgrp$lb.term)
  bndgrp$ub.term <- ifelse(bndgrp$ub == COBRAR_SETTINGS("MAXIMUM"),"default_ub", bndgrp$ub.term)
  bndgrp_para <- data.frame(bnd = c(bndgrp$lb.term, bndgrp$ub.term),
                            val = c(bndgrp$lb, bndgrp$ub))
  bndgrp_para <- bndgrp_para[!duplicated(bndgrp_para$bnd),]
  bndgrp_para$SBO <- ifelse(grepl("^default_",bndgrp_para$bnd),626,625)
  bndgrp_para <- bndgrp_para[order(bndgrp_para$bnd),]

  out <- writeSBML(
    file_path = file_path,

    # Model fields
    mod_id = model@mod_id,
    mod_name = model@mod_name,
    mod_desc = model@mod_desc,

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
    react_rev = ifelse(bndgrp$lb < 0 & bndgrp$ub > 0, TRUE, FALSE)
  )

  return(out)
}
