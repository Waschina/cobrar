#' FROG: Reproducible fitness and robustness diagnostics for constraint-based
#' models
#'
#' Run a standardized FROG analysis on a model including (1) Objective Function
#' Values, (2) Flux Variability Analysis (FVA), (3) Gene Deletion Fluxes, (4)
#' Reaction Deletion Fluxes.
#'
#' @param model Model of class \link{ModelOrg} or path to SBML file.
#' @param outdir Path to output directory, where report files will be saved.
#' @param ... Additional parameters passed on to the specific method instance.
#' 
#' @returns TRUE if frog report export was successful.
#'
#' @examples
#' fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
#' mod <- readSBMLmod(fpath)
#' 
#' output_dir <- tempdir()
#' frog(mod, outdir = output_dir)
#' 
#' # list files in output directory
#' dir(output_dir)
#' 
#' @references
#' K. Raman et al., “FROG Analysis Ensures the Reproducibility of Genome Scale
#' Metabolic Models.,” bioRxiv, Sept. 2024. doi: 10.1101/2024.09.24.614797.
#' 
#' @docType methods
#' @rdname frog-methods
#' @export
setGeneric("frog", valueClass = "logical", function(model, ...) {
  standardGeneric("frog")
})

#' @rdname frog-methods
#' @aliases frog,ModelOrg
setMethod("frog", signature(model = "ModelOrg"),
          function(model, outdir = "frog_report") {
            
            # Check if output directory is valid
            if(!dir.exists(outdir) ) {
              dir_success <- dir.create(outdir)
              if(!dir_success) {
                stop(paste0("Output directory '", outdir,
                            "' cannot be created."))
              }
            }
            
            xmlfile <- paste0(outdir,"/", model@mod_id, ".xml")
            
            # run frog
            out_frog <- run_frog(model, filename = basename(xmlfile))
            
            # export results and model
            xml_success <- writeSBMLmod(model, xmlfile)
            export_success <- export_frog_report(out_frog, outdir, xmlfile)

            if(!xml_success)
              warning("Exporting model as SBML failed.")
            
            if(!export_success)
              warning("Writing FROG report failed.")
            
            return(all(c(xml_success, export_success)))
          }
)

#' @rdname frog-methods
#' @aliases frog,character
setMethod("frog", signature(model = "character"),
          function(model, outdir = "frog_report") {
            
            # Check if output directory is valid
            if(!dir.exists(outdir) ) {
              dir_success <- dir.create(outdir)
              if(!dir_success) {
                stop(paste0("Output directory '", outdir,
                            "' cannot be created."))
              }
            }
            
            xmlfile <- model
            model <- readSBMLmod(xmlfile)
            
            # run frog
            out_frog <- run_frog(model, filename = basename(xmlfile))
            
            # export results
            export_success <- export_frog_report(out_frog, outdir, xmlfile)
            
            if(!export_success)
              warning("Writing FROG report failed.")
            
            return(export_success)
          }
)

run_frog <- function(model, filename) {
  
  out_frog <- list()
  
  # (1) FBA of objective function
  sol <- fba(model)
  out_frog$fba <- data.frame(model = filename,
                             objective = printObjFunc(model),
                             status = sol@stat_term,
                             value = sol@obj)
  
  # (2) FVA
  sol_fva <- fva(model)
  out_frog$fva <- data.frame(
    model = filename,
    objective = printObjFunc(model),
    reaction = sol_fva$react,
    flux = sol@fluxes,
    status = sol@stat_term,
    minimum = sol_fva$min.flux,
    maximum = sol_fva$max.flux
  )
  
  # (3) Gene Deletion Fluxes
  
  # create warm start object
  LPprob <- new(paste0("LPproblem_",COBRAR_SETTINGS("SOLVER")),
                name = paste0("LP_", model@mod_id),
                method = COBRAR_SETTINGS("METHOD"))
  
  loadLPprob(LPprob,
             nCols = react_num(model),
             nRows = met_num(model)+constraint_num(model),
             mat   = rbind(model@S, model@constraints@coeff),
             ub    = ifelse(abs(model@uppbnd)>COBRAR_SETTINGS("MAXIMUM"),
                            sign(model@uppbnd)*COBRAR_SETTINGS("MAXIMUM"),
                            model@uppbnd),
             lb    = ifelse(abs(model@lowbnd)>COBRAR_SETTINGS("MAXIMUM"),
                            sign(model@lowbnd)*COBRAR_SETTINGS("MAXIMUM"),
                            model@lowbnd),
             obj   = model@obj_coef,
             rlb   = c(rep(0, met_num(model)),
                       model@constraints@lb),
             rtype = c(rep("E", met_num(model)),
                       model@constraints@rtype),
             lpdir = substr(model@obj_dir,1,3),
             rub   = c(rep(NA, met_num(model)),
                       model@constraints@ub),
             ctype = NULL
  )
  
  if(gene_num(model) > 0) {
    res3 <- numeric(gene_num(model))
    res3_stat <- character(gene_num(model))
    for(j in 1:gene_num(model)) {
      i <- react_pos(model, geneDel(model, j))
      if(length(i) > 0) {
        bu_lp <- model@lowbnd[i]
        bu_up <- model@uppbnd[i]
        bu_objc <- model@obj_coef[i]
        n <- length(i)
        
        setColsBndsObjCoefs(LPprob, i, lb = rep(0,n), ub = rep(0,n), obj_coef = bu_objc)
        lp_ok   <- solveLp(LPprob)
        lp_stat <- getSolStat(LPprob)
        
        res3[j] <- getObjValue(LPprob)
        res3_stat[j] <- lp_stat$term
        
        setColsBndsObjCoefs(LPprob, i, lb = bu_lp, ub = bu_up, obj_coef = bu_objc)
        
      } else {
        res3[j] <- sol@obj
        res3_stat[j] <- sol@stat_term
      }
    }
    out_frog$gene_del <- data.frame(
      model = filename,
      objective = printObjFunc(model),
      gene = model@allGenes,
      status = res3_stat,
      value = res3
    )
  } else {
    out_frog$gene_del <- data.frame(
      model = character(0L),
      objective = character(0L),
      gene = character(0L),
      status = character(0L),
      value = numeric(0L)
    )
  }
  
  # (4) Reaction deletion fluxes
  
  res4 <- numeric(react_num(model))
  res4_stat <- character(react_num(model))
  for(i in 1:react_num(model)) {
    bu_lp <- model@lowbnd[i]
    bu_up <- model@uppbnd[i]
    bu_objc <- model@obj_coef[i]
    
    setColsBndsObjCoefs(LPprob, i, lb = 0, ub = 0, obj_coef = bu_objc)
    lp_ok   <- solveLp(LPprob)
    lp_stat <- getSolStat(LPprob)
  
    res4[i] <- getObjValue(LPprob)
    res4_stat[i] <- lp_stat$term
    
    setColsBndsObjCoefs(LPprob, i, lb = bu_lp, ub = bu_up, obj_coef = bu_objc)
  }
  
  out_frog$react_del <- data.frame(
    model = filename,
    objective = printObjFunc(model),
    reaction = model@react_id,
    status = res4_stat,
    value = res4
  )
  
  out_frog
}

#' @importFrom jsonlite write_json
#' @importFrom utils write.table packageVersion
export_frog_report <- function(frogres, outdir, xmlfile) {
  
  export_success <- c()
  
  # fluxes
  export_success[1] <- write.table(frogres$fba,
                                   paste0(outdir, "/01_objective.tsv"),
                                   sep = "\t", quote = FALSE, row.names = FALSE)
  export_success[2] <- write.table(frogres$fva,
                                   paste0(outdir, "/02_fva.tsv"),
                                   sep = "\t", quote = FALSE, row.names = FALSE)
  export_success[3] <- write.table(frogres$gene_del,
                                   paste0(outdir,"/03_gene_deletion.tsv"),
                                   sep = "\t", quote = FALSE, row.names = FALSE)
  export_success[4] <- write.table(frogres$react_del,
                                   paste0(outdir,"/03_reaction_deletion.tsv"),
                                   sep = "\t", quote = FALSE, row.names = FALSE)
  
  # meta
  metafrog <- list(
    software = list(
      frog = list(
        name = "cobrar FROG",
        version = as.character(packageVersion("cobrar")),
        url = "https://github.com/Waschina/cobrar/R/frog.R"
      ),
      toolbox = list(
        name = "cobrar",
        version = as.character(packageVersion("cobrar")),
        url = "https://github.com/Waschina/cobrar"
      ),
      solver = list(
        name = COBRAR_SETTINGS("SOLVER"),
        version = getGLPKVersion()
      )
    ),
    model_filename = basename(xmlfile),
    frog_date = Sys.Date(),
    # model_sha256 = sha256sum(xmlfile), # not available in 'tools' < 4.5.1
    model_md5 = md5sum(xmlfile),
    environment = R.version$version.string,
    frog_version = "0.1.4"
  )
  
  export_success[5] <-write_json(metafrog, 
                                 paste0(outdir,"/metadata.json"),
                                 auto_unbox = TRUE, pretty = TRUE)
  
  all(export_success)
}
