test_that("nbcomp.bootspls returns CV matrix and sensible optima", {
  skip_on_cran()
  skip_if_not_installed("bootPLS")
  set.seed(20240912)
  # Small synthetic regression
  X <- matrix(rnorm(40*6), nrow = 40, ncol = 6)
  beta <- c(1, -0.5, 0.8, rep(0,3))
  y <- as.numeric(X %*% beta + rnorm(40, sd = 0.5))
  eta <- c(0.2, 0.6)
  res <- nbcomp.bootspls(x = X, y = y, eta = eta, R = 60, maxnt = 2,
                         plot.it = FALSE, typeBCa = FALSE, verbose = FALSE)
  expect_type(res, "list")
  expect_true(all(c("mspemat","eta.opt","K.opt") %in% names(res)))
  expect_equal(nrow(res$mspemat), length(eta))
  expect_true(res$eta.opt %in% eta)
  expect_true(res$K.opt %in% 1:2)
  expect_true(all(is.finite(res$mspemat)))
})
