test_that("Package loads quietly", {
  expect_true(requireNamespace("bootPLS", quietly = TRUE))
})
