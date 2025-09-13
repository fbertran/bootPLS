test_that("signpred2 can draw without error on a small binary matrix", {
  skip_on_cran()
  skip_if_not_installed("bootPLS")
  set.seed(42)
  simbin <- matrix(rbinom(30, 1, .3), nrow = 12, ncol = 5)
  colnames(simbin) <- paste0("C", 1:ncol(simbin))
  rownames(simbin) <- paste0("R", 1:nrow(simbin))
  # just ensure it does not error; plot to null device
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(suppressWarnings(signpred2(simbin, labsize = .8, plotsize = 4)))
})
