test_that("deadEndMetabolites identifies blocked metabolites", {
  model <- make_dead_end_model()

  result <- deadEndMetabolites(model)

  expect_named(result, c("dem", "der"))
  expect_equal(result$dem, "c")
  expect_equal(result$der, "a_to_c")
})
