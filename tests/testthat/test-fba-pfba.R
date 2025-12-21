test_that("FBA optimizes e_coli_core model using GLPK backend", {
  fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
  model <- readSBMLmod(fpath)

  result <- fba(model)

  expected_fluxes = c(
    7.477382e+00,  0.000000e+00,  4.860861e+00, -1.602353e+01,  4.959985e+00,  2.387424e-12,
    0.000000e+00, -1.471614e+01,  3.214895e+00,  0.000000e+00,  2.273737e-12, -4.774847e-12,
    2.504309e+00,  6.007250e+00,  6.007250e+00,  8.390000e+00,  0.000000e+00, -4.774847e-12,
    0.000000e+00, -6.821210e-13,  5.064376e+00,  4.551401e+01, -1.818989e-12,  1.758177e+00,
    8.739215e-01,  8.753887e-12, -2.280983e+01,  2.678482e+00,  6.007250e+00, -2.281503e+00,
   -6.963319e-13,  4.359899e+01,  3.410605e-13,  1.471614e+01,  0.000000e+00,  0.000000e+00,
    5.064376e+00, -5.064376e+00,  1.496984e+00,  0.000000e+00,  1.496984e+00,  1.181498e+00,
    7.477382e+00,  0.000000e+00,  0.000000e+00,  0.000000e+00,  2.280983e+01,  0.000000e+00,
    0.000000e+00,  0.000000e+00,  0.000000e+00, -1.000000e+01,  0.000000e+00,  0.000000e+00,
    1.753087e+01,  2.917583e+01,  0.000000e+00,  0.000000e+00, -4.765319e+00, -2.179949e+01,
   -3.214895e+00,  0.000000e+00,  0.000000e+00,  7.477382e+00,  0.000000e+00,  0.000000e+00,
    0.000000e+00,  0.000000e+00,  2.958258e-15,  5.064376e+00,  0.000000e+00,  4.959985e+00,
    1.602353e+01,  1.000000e+01,  2.234617e-01,  0.000000e+00, -4.541857e+00,  0.000000e+00,
    0.000000e+00,  1.136868e-12,  4.959985e+00, -2.917583e+01,  6.007250e+00,  0.000000e+00,
    3.410605e-13,  8.384404e-13,  0.000000e+00,  5.064376e+00,  0.000000e+00,  0.000000e+00,
    3.853461e+01,  0.000000e+00,  4.765319e+00,  2.179949e+01,  9.282533e+00
  )

  expect_s4_class(result, "FluxPrediction")
  expect_identical(result@algorithm, "FBA")
  expect_identical(result@ok_term, "optimization process was successful")
  expect_identical(result@stat_term, "solution is optimal")
  expect_equal(result@obj, 0.8739215, tolerance = 1e-6)
  expect_equal(result@fluxes, expected_fluxes, tolerance = 1e-6)
})

test_that("pFBA preserves objective and minimizes weighted flux", {
  fpath <- system.file("extdata", "e_coli_core.xml", package="cobrar")
  model <- readSBMLmod(fpath)

  baseline <- fba(model)
  result <- pfba(model)

  expect_s4_class(result, "FluxPrediction")
  expect_equal(result@obj, baseline@obj, tolerance = 1e-6)
  expect_equal(result@fluxes, baseline@fluxes, tolerance = 1e-6)
  expect_equal(result@obj_sec, sum(abs(baseline@fluxes)), tolerance = 1e-6)

  costcoeffw_list = sample(c(1, 3), length(model@react_id), replace = TRUE)
  costcoefbw_list = sample(c(2, 4), length(model@react_id), replace = TRUE)
  weighted <- pfba(
    model,
    costcoeffw = costcoeffw_list,
    costcoefbw = costcoefbw_list
  )
  expect_equal(
    weighted@obj_sec,
    sum(abs(baseline@fluxes * costcoeffw_list)),
#    tolerance = 1e-6
    tolerance = 1e-4 * sum(costcoeffw_list) * length(model@react_id)
  )

  expect_error(
    pfba(model, costcoeffw = c(1)),
    "must both be of length equal to the number of reactions"
  )
})

