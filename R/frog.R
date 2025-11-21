# FROG analysis for CobraR
# Public API (do not break): frog_run(), frog_fva(), frog_reaction_deletion(), frog_gene_deletion(),
# frog_objective_checks(), frog_export()
# Depends: cobrar (>= 0.2.2), jsonlite, utils, tools
#' FROG: Reproducible fitness & robustness diagnostics for constraint-based models
#'
#' Run a standardized FROG analysis on a CobraR model and return a structured list
#' containing FBA baseline, FVA intervals, reaction and gene single-deletion scans,
#' objective sanity checks, and metadata. Pair with [frog_export()] to write a
#' portable FROG report (CSVs + JSON; optional COMBINE archive `.omex`).
#'
#' @param mod A CobraR model (ModelOrg) loaded via `cobrar::readSBMLmod()`.
#' @param obj_rxn Optional reaction id to set as objective (default: keep current).
#' @param sense "max" or "min" objective direction (default: "max").
#' @param fva_frac Fraction of optimal objective for FVA (default: 1.0).
#' @param ko_threads Integer number of worker processes for deletion scans
#'   (use >1 on Unix with `parallel`; Windows runs serial).
#' @return A list of class `frog_result` with elements: `fba`, `fva`,
#'   `reaction_deletion`, `gene_deletion`, `objective_checks`, `meta`.
#' @examples
#' \dontrun{
#' f <- system.file("extdata","e_coli_core.xml", package="cobrar")
#' mod <- cobrar::readSBMLmod(f)
#' res <- frog_run(mod, obj_rxn = "BIOMASS_Ecoli_core_w_GAM", sense = "max",
#'                 fva_frac = 1.0, ko_threads = 1)
#' frog_export(res, out_dir = "FROG_report", make_omex = TRUE)
#' }
#' @export

frog_run <- function(mod, obj_rxn = NULL, sense = "max",
                     fva_frac = 1.0, ko_threads = 1) {
  if (!requireNamespace("cobrar", quietly = TRUE)) {
    stop("Package 'cobrar' is required.")
  }
  
  sense <- match.arg(sense)
  
  # set objective if requested (non-destructive: work on a copy)
  if (!is.null(obj_rxn)) {
    rid <- tryCatch(cobrar::react_pos(mod, obj_rxn), error=function(e) integer())
    if (length(rid) != 1L) stop("Objective reaction not found: ", obj_rxn)
    obj <- rep(0, tryCatch(cobrar::react_num(mod), error=function(e) 0L))
    obj[rid] <- 1
    mod <- cobrar::changeObjFunc(mod, obj_rxn)
  }
  mod <- cobrar::setObjDir(mod, ifelse(sense=="max", "maximize", "minimize"))

  # FBA reference
  fba_res <- suppressMessages(cobrar::fba(mod))
  if (is.na(fba_res@obj)) stop("FBA failed: no optimal solution")

  # FVA
  fva_res <- frog_fva(mod, fba_res, fva_frac)

  # Reaction deletion
  rxn_del <- frog_reaction_deletion(mod, ko_threads = ko_threads)

  # Gene deletion (robust fallback if mapping isn't available)
  gene_del <- frog_gene_deletion(mod, ko_threads = ko_threads)

  # Objective checks
  obj_checks <- frog_objective_checks(mod, fba_res)

  # Meta/env
  meta <- frog_metadata(mod)

  out <- list(
    fba = frog_tidy_fba(fba_res),
    fva = fva_res,
    reaction_deletion = rxn_del,
    gene_deletion = gene_del,
    objective_checks = obj_checks,
    meta = meta
  )
  class(out) <- c("frog_result","list")
  out
}

#' Run Flux Variability Analysis for FROG
#' @param mod CobraR model
#' @param fba_res Optional FBA result to reuse
#' @param fva_frac Fraction of optimal objective (0,1]
#' @return data.frame with columns: reaction, min, max
#' @export
frog_fva <- function(mod, fba_res = NULL, fva_frac = 1.0) {
  if (!is.numeric(fva_frac) || fva_frac <= 0 || fva_frac > 1) 
    stop("fva_frac must be in (0,1].")
    if (is.null(fba_res)) fba_res <- cobrar::fba(mod)
  if (is.na(fba_res@obj)) stop("No optimal solution for FVA baseline")
  suppressMessages({
    fv <- cobrar::fva(mod, opt.factor = fva_frac)})
  if(!all(mod@react_id == fv$react)) stop ('Reaction ids in fva dont match reaction ids in fba')
  data.frame(
    model =   mod@mod_id,
    objective = 'obj',
    reaction = fv$react,
    flux = fba_res@fluxes,
    status = fba_res@stat_term,
    minimum = fv$min,
    maximum = fv$max,
    fraction_optimum = fva_frac,
    stringsAsFactors = FALSE
  )
}

frog_tidy_fba <- function(fba_res) {
  data.frame(
    model = mod@mod_id,
    objective = 'obj',
    # reaction = names(fba_res@fluxes),
    # flux = as.numeric(fba_res@fluxes),
    solver_status = fba_res@stat_term,
    value = fba_res@obj,
    stringsAsFactors = FALSE
  )
}

#' Reaction single-deletion scan
#' @param mod CobraR model
#' @param ko_threads parallel workers (Unix only)
#' @return data.frame: reaction, objective, lethal
#' @export
frog_reaction_deletion <- function(mod, ko_threads = 1) {
  #rxns <- tryCatch(cobrar::react_pos(mod), error=function(e) integer())
  rxns <- mod@react_id
  #base_obj <- tryCatch(cobrar::fba(mod)@obj, error=function(e) NA_real_)
  base_obj <- suppressMessages(cobrar::fba(mod))@obj
  worker <- function(i) {
    m <- mod
    m <- cobrar::rmReact(m, i)
    fr <- suppressMessages(cobrar::fba(m))
    data.frame(
      reaction = tryCatch(cobrar::react_pos(mod, i, to = "id"), error=function(e) as.character(i)),
      objective = ifelse(is.na(fr@obj), NA_real_, fr@obj),
      solver_status = fr@stat_term,
      lethal = ifelse(is.na(fr@obj) || is.na(base_obj), NA, fr@obj < 1e-9 & base_obj > 1e-9),
      stringsAsFactors = FALSE
    )
  }
  if (length(rxns) == 0L) return(data.frame(reaction=character(), objective=numeric(), solver_status =character(), lethal=logical()))
  if (.Platform$OS.type != "windows" && ko_threads > 1 && requireNamespace("parallel", quietly = TRUE)) {
    parts <- parallel::mclapply(rxns, worker, mc.cores = ko_threads)
  } else {
    parts <- lapply(rxns, worker)
  }
  do.call(rbind, parts)
}



#' @param mod CobraR model
#' @param ko_threads parallel workers (Unix only)
#' @return data.frame: gene, objective, lethal
#' @export

####Do we want to include unique genes? (below code from MATLAB)
# if (uniqueGene == 1)
#   
#   % detect whether there are alternate transcripts
# transcriptsPresent = false;
# if any(~cellfun(@isempty, regexp(model.genes,'\.[0-9]+$'))) %If there are any genes that end on a transcript.
# transcriptsPresent = true;
# [geneList,rem] = strtok(geneList,'.');
# geneList = unique(geneList);
# nGenes = length(geneList);
# nDelGenes = length(geneList);
# else
#   nGenes = length(model.genes);
# nDelGenes = length(geneList);

frog_gene_deletion <- function(mod, ko_threads = 1) {
  g_ids <- mod@allGenes
 
  base_obj <- tryCatch(cobrar::fba(mod)@obj, error=function(e) NA_real_)
  # Define worker function using the gene_reaction_map
  worker <- function(g) {
    m <- mod
     m <- tryCatch(
        cobrar::rmGene(m, g),
        error = function(e) m
      )
    
    # Run FBA
    fr <- tryCatch(suppressMessages(cobrar::fba(m)), error = function(e) NULL)
    
    # Return result row
    data.frame(
      gene = g,
      objective = if (!is.null(fr)) fr@obj else NA_real_,
      solver_status = fr@stat_term,
      lethal = if (is.null(fr) || is.na(base_obj)) NA else (fr@obj < 1e-9 & base_obj > 1e-9),
      stringsAsFactors = FALSE
    )
  }
  if (length(g_ids) == 0L) {
    return(data.frame(gene=character(), objective=numeric(), solver_status = character() ,lethal=logical()))
  }
  if (.Platform$OS.type != "windows" && ko_threads > 1 && requireNamespace("parallel", quietly = TRUE)) {
    parts <- parallel::mclapply(g_ids, worker, mc.cores = ko_threads)
  } else {
    parts <- lapply(g_ids, worker)
  }
  do.call(rbind, parts)
}

#' Objective sanity checks
#' @param mod CobraR model
#' @param fba_res optional FBA result
#' @return named list with feasibility, objective, direction, etc.
#' @export
frog_objective_checks <- function(mod, fba_res = NULL) {
  if (is.null(fba_res)) fba_res <- cobrar::fba(mod)
  obj <- fba_res@obj
  list(
    feasible = !is.na(obj),
    objective_value = obj,
    objective_direction = tryCatch(mod@obj_dir, error=function(e) NA),
    nonzero_flux_at_objective = any(abs(fba_res@fluxes) > 1e-12),
    biomass_guess = tryCatch({
      bm <- cobrar::guessBMReaction(mod)
      if (length(bm) > 0) bm else NA_character_
    }, error = function(e) NA_character_)
  )
}

frog_metadata <- function(mod) {
  list(
    timestamp_utc = as.character(Sys.time()),
    cobrar_version = tryCatch(as.character(utils::packageVersion("cobrar")), error=function(e) NA),
    sbml_level = 3, sbml_version = 2, fbc_version = 2,
    model = list(
      id = tryCatch(mod@mod_id, error=function(e) NA),
      name = tryCatch(mod@mod_name, error=function(e) NA),
      reactions = tryCatch(cobrar::react_num(mod), error=function(e) NA),
      metabolites = tryCatch(cobrar::met_num(mod), error=function(e) NA),
      genes = tryCatch(cobrar::gene_num(mod), error=function(e) NA)
    ),
    checksums = list()
  )
}

#' Export FROG results to a portable report
#'
#' Writes CSVs (fba, fva, reaction_deletion, gene_deletion, objective_checks),
#' a JSON summary (frog.json), and optionally packages them as `FROG.omex`
#' (a COMBINE archive with a minimal `omex-manifest.xml`).
#'
#' @param frog A `frog_result` list from [frog_run()].
#' @param out_dir Output directory (created if missing).
#' @param make_omex Logical, also create `FROG.omex` one level above `out_dir`.
#' @return Invisibly returns path to OMEX or directory.
#' @export
frog_export <- function(frog, out_dir = "FROG_report", make_omex = TRUE) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Package 'jsonlite' required for export.")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # CSVs
  utils::write.csv(frog$fba, file.path(out_dir, "fba.csv"), row.names = FALSE)
  utils::write.csv(frog$fva, file.path(out_dir, "fva.csv"), row.names = FALSE)
  utils::write.csv(frog$reaction_deletion, file.path(out_dir, "reaction_deletion.csv"), row.names = FALSE)
  utils::write.csv(frog$gene_deletion, file.path(out_dir, "gene_deletion.csv"), row.names = FALSE)
  utils::write.csv(
    data.frame(key = names(frog$objective_checks), value = I(unlist(frog$objective_checks))),
    file.path(out_dir, "objective_checks.csv"), row.names = FALSE
  )

  # JSON
  meta <- frog$meta
  meta$checksums <- list(
    fba_csv = tryCatch(tools::md5sum(file.path(out_dir, "fba.csv"))[[1]], error=function(e) NA),
    fva_csv = tryCatch(tools::md5sum(file.path(out_dir, "fva.csv"))[[1]], error=function(e) NA),
    reaction_deletion_csv = tryCatch(tools::md5sum(file.path(out_dir, "reaction_deletion.csv"))[[1]], error=function(e) NA),
    gene_deletion_csv = tryCatch(tools::md5sum(file.path(out_dir, "gene_deletion.csv"))[[1]], error=function(e) NA),
    objective_checks_csv = tryCatch(tools::md5sum(file.path(out_dir, "objective_checks.csv"))[[1]], error=function(e) NA)
  )
  jsonlite::write_json(
    list(schema = "https://biomodels.net/frog/schema/0.1", meta = meta),
    file.path(out_dir, "frog.json"), pretty = TRUE, auto_unbox = TRUE
  )

  # OMEX manifest
  manifest <- paste0(
    '<?xml version="1.0" encoding="utf-8"?>\n',
    '<omexManifest xmlns="http://identifiers.org/combine.specifications/omex-manifest">\n',
    '  <content location="./" format="http://identifiers.org/combine.specifications/omex"/>\n',
    '  <content location="frog.json" format="application/json"/>\n',
    '  <content location="fba.csv" format="text/csv"/>\n',
    '  <content location="fva.csv" format="text/csv"/>\n',
    '  <content location="reaction_deletion.csv" format="text/csv"/>\n',
    '  <content location="gene_deletion.csv" format="text/csv"/>\n',
    '  <content location="objective_checks.csv" format="text/csv"/>\n',
    '</omexManifest>\n'
  )
  writeLines(manifest, file.path(out_dir, "omex-manifest.xml"))
  if (isTRUE(make_omex)) {
    omex_name <- "FROG.omex"
    omex_tmp_path <- file.path(tempdir(), omex_name)
    omex_final_path <- normalizePath(file.path(out_dir, omex_name), mustWork = FALSE)
    
    old <- getwd(); on.exit(setwd(old), add = TRUE)
    setwd(out_dir)
    
    # Create OMEX archive inside tempdir
    utils::zip(zipfile = omex_tmp_path,
               files = c("frog.json","fba.csv","fva.csv",
                         "reaction_deletion.csv","gene_deletion.csv",
                         "objective_checks.csv","omex-manifest.xml"))
    
    # Restore working directory
    setwd(old)
    
    # Clean report folder
    unlink(list.files(out_dir, full.names = TRUE))
    
    # Move OMEX into report folder
    file.copy(omex_tmp_path, omex_final_path, overwrite = TRUE)
    unlink(omex_tmp_path)
    
    message("Wrote ", omex_final_path)
    return(invisible(omex_final_path))
  }
  invisible(out_dir)

}
