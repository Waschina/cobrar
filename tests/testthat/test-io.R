test_that("SBML import loads expected metadata", {
  sbml_path <- system.file("extdata", "e_coli_core.xml", package = "cobrar")
  skip_if(identical(sbml_path, ""), "SBML fixture not installed")

  model <- readSBMLmod(sbml_path)

  expect_s4_class(model, "ModelOrg")
  expect_equal(model@mod_id, "e_coli_core")
  expect_equal(model@obj_dir, "maximize")
  expect_true("BIOMASS_Ecoli_core_w_GAM" %in% model@react_id)
  expect_equal(
    model@obj_coef[match("BIOMASS_Ecoli_core_w_GAM", model@react_id)],
    1
  )
  expect_equal(length(model@mod_compart), 2)
})

test_that("SBML export/import round-trips toy model", {
  toy <- make_minimal_flux_model()
  tmp <- tempfile(fileext = ".xml")
  on.exit(unlink(tmp), add = TRUE)

  expect_true(writeSBMLmod(toy, tmp))
  expect_true(file.exists(tmp))
  
  sol <- fba(toy)

  roundtrip <- readSBMLmod(tmp)

  expect_equal(roundtrip@mod_id,   toy@mod_id)
  expect_equal(roundtrip@obj_dir,  toy@obj_dir)
  expect_equal(as.matrix(roundtrip@S), as.matrix(toy@S))
  expect_equal(roundtrip@lowbnd,   toy@lowbnd)
  expect_equal(roundtrip@uppbnd,   toy@uppbnd)
  expect_equal(roundtrip@obj_coef, toy@obj_coef)
  expect_equal(roundtrip@react_id, toy@react_id)
  expect_equal(roundtrip@met_id,   toy@met_id)
  
  rt_sol <- fba(roundtrip)
  expect_equal(rt_sol@obj, sol@obj)
  expect_identical(rt_sol@ok, sol@ok)
})
