# tests for FROG functions
testthat::test_that("FROG pipeline runs on E. coli core", {
  testthat::skip_if_not_installed("cobrar")
  f <- system.file("extdata","e_coli_core.xml", package="cobrar")
  testthat::skip_if(!file.exists(f), "e_coli_core.xml not available")
  mod <- cobrar::readSBMLmod(f)

  res <- frog_run(mod,
                  obj_rxn  = "BIOMASS_Ecoli_core_w_GAM",
                  sense    = "max",
                  fva_frac = 1.0,
                  ko_threads = 1)

  testthat::expect_s3_class(res, "frog_result")
  testthat::expect_true(is.data.frame(res$fba))
  testthat::expect_true(is.data.frame(res$fva))
  testthat::expect_true(is.data.frame(res$reaction_deletion))
  testthat::expect_true(is.data.frame(res$gene_deletion))
  testthat::expect_true(is.list(res$objective_checks))

  # Objective should be numeric; allow zero/nonzero depending on constraints
  obj <- unique(res$fba$objective)
  testthat::expect_true(is.numeric(obj))
})
