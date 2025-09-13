test_that("nbcomp.bootsgpls works on small dataset and returns expected fields", {
  skip_on_cran()
  skip_if_not_installed("spls")
  data(prostate, package = "spls")
  set.seed(20240912)
  # keep it tiny for speed
  X <- as.matrix((prostate$x)[, 1:20])
  y <- as.numeric(prostate$y)
  eta <- 0.2
  res <- nbcomp.bootsgpls(x = X, y = y, R = 50, eta = eta, maxnt = 2,
                          plot.it = FALSE, typeBCa = FALSE, verbose = FALSE)
  expect_type(res, "list")
  expect_true(all(c("err.mat","eta.opt","K.opt") %in% names(res)))
  expect_equal(nrow(res$err.mat), 1L)
  expect_true(res$eta.opt == eta)
  expect_true(res$K.opt %in% 1:2)
})
