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
frog_run <- function(mod, obj_rxn = NULL, sense = c("max","min"),
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
    mod <- cobrar::changeObjFunc(mod, obj)
  }
  mod <- cobrar::setObjDir(mod, ifelse(sense=="max", "max", "min"))

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
#' @param fraction_of_opt Fraction of optimal objective (0,1]
#' @return data.frame with columns: reaction, min, max
#' @export
frog_fva <- function(mod, fba_res = NULL, fraction_of_opt = 1.0) {
  if (!is.numeric(fraction_of_opt) || fraction_of_opt <= 0 || fraction_of_opt > 1)
    stop("fraction_of_opt must be in (0,1].")
  if (is.null(fba_res)) fba_res <- cobrar::fba(mod)
  if (is.na(fba_res@obj)) stop("No optimal solution for FVA baseline")
  suppressMessages({
    fv <- cobrar::fva(mod, fraction_of_optimum = fraction_of_opt)
  })
  data.frame(
    reaction = fv$reaction,
    min = fv$min,
    max = fv$max,
    stringsAsFactors = FALSE
  )
}

frog_tidy_fba <- function(fba_res) {
  data.frame(
    reaction = names(fba_res@fluxes),
    flux = as.numeric(fba_res@fluxes),
    objective = fba_res@obj,
    solver_status = fba_res@stat,
    stringsAsFactors = FALSE
  )
}

#' Reaction single-deletion scan
#' @param mod CobraR model
#' @param ko_threads parallel workers (Unix only)
#' @return data.frame: reaction, objective, lethal
#' @export
frog_reaction_deletion <- function(mod, ko_threads = 1) {
  rxns <- tryCatch(cobrar::react_pos(mod), error=function(e) integer())
  base_obj <- tryCatch(cobrar::fba(mod)@obj, error=function(e) NA_real_)
  worker <- function(i) {
    m <- mod
    # KO by fixing bounds to 0
    m <- cobrar::changeBounds(m, i, lb = 0, ub = 0)
    fr <- suppressMessages(cobrar::fba(m))
    data.frame(
      reaction = tryCatch(cobrar::react_pos(mod, i, to = "id"), error=function(e) as.character(i)),
      objective = ifelse(is.na(fr@obj), NA_real_, fr@obj),
      lethal = ifelse(is.na(fr@obj) || is.na(base_obj), NA, fr@obj < 1e-9 & base_obj > 1e-9),
      stringsAsFactors = FALSE
    )
  }
  if (length(rxns) == 0L) return(data.frame(reaction=character(), objective=numeric(), lethal=logical()))
  if (.Platform$OS.type != "windows" && ko_threads > 1 && requireNamespace("parallel", quietly = TRUE)) {
    parts <- parallel::mclapply(rxns, worker, mc.cores = ko_threads)
  } else {
    parts <- lapply(rxns, worker)
  }
  do.call(rbind, parts)
}

#' Gene single-deletion scan (via GPR where available)
#' @param mod CobraR model
#' @param ko_threads parallel workers (Unix only)
#' @return data.frame: gene, objective, lethal
#' @export
frog_gene_deletion <- function(mod, ko_threads = 1) {
  g_ids <- tryCatch(cobrar::gene_pos(mod, to = "id"), error=function(e) character())
  # Preferred mapping via cobrar::geneDel() if available
  gmap <- tryCatch(cobrar::geneDel(mod), error=function(e) NULL)

  base_obj <- tryCatch(cobrar::fba(mod)@obj, error=function(e) NA_real_)
  worker <- function(g) {
    m <- mod
    affected <- character()
    if (!is.null(gmap)) {
      # expect data.frame with columns $gene and $reaction; be permissive
      gg <- tryCatch(gmap$gene, error=function(e) NULL)
      rr <- tryCatch(gmap$reaction, error=function(e) NULL)
      if (!is.null(gg) && !is.null(rr)) {
        affected <- unique(rr[gg == g])
      }
    }
    if (length(affected) > 0) {
      m <- cobrar::changeBounds(m, cobrar::react_pos(m, affected), lb = 0, ub = 0)
    }
    fr <- suppressMessages(cobrar::fba(m))
    data.frame(
      gene = g,
      objective = ifelse(is.na(fr@obj), NA_real_, fr@obj),
      lethal = ifelse(is.na(fr@obj) || is.na(base_obj), NA, fr@obj < 1e-9 & base_obj > 1e-9),
      stringsAsFactors = FALSE
    )
  }
  if (length(g_ids) == 0L) {
    return(data.frame(gene=character(), objective=numeric(), lethal=logical()))
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
    objective_direction = tryCatch(cobrar::COBRAR_SETTINGS()$obj_dir, error=function(e) NA),
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
      id = tryCatch(attr(mod, "id"), error=function(e) NA),
      name = tryCatch(attr(mod, "name"), error=function(e) NA),
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
    omex <- file.path(normalizePath(file.path(out_dir, "..")), "FROG.omex")
    old <- getwd(); on.exit(setwd(old), add = TRUE)
    setwd(out_dir)
    utils::zip(zipfile = omex, files = c("frog.json","fba.csv","fva.csv",
                                         "reaction_deletion.csv","gene_deletion.csv",
                                         "objective_checks.csv","omex-manifest.xml"))
    message("Wrote ", omex)
    return(invisible(omex))
  }
  invisible(out_dir)
}
