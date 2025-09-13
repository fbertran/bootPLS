test_that("nbcomp.bootplsR selects a valid number of components", {
  skip_on_cran()
  skip_if_not_installed("plsRglm")
  set.seed(20240912)
  data(pine, package = "plsRglm")
  X <- as.matrix(pine[, 1:8])
  y <- as.numeric(log(pine[, 11]))
  k <- nbcomp.bootplsR(Y = y, X = X, R = 60, sim = "ordinary",
                       ncpus = 1, parallel = "no", typeBCa = FALSE, verbose = FALSE)
  expect_true(is.numeric(k))
  expect_true(abs(k - round(k)) < .Machine$double.eps^0.5) # integer-ish
  expect_true(k >= 0 && k <= 8)
})
