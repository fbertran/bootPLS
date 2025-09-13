test_that("simul_data_UniYX_gamma generates expected named vector and can be replicated", {
  skip_on_cran()
  skip_if_not_installed("bootPLS")
  set.seed(314)
  # ensure small and deterministic parameters
  res1 <- simul_data_UniYX_gamma(totdim = 20, ncomp = 2, lvar = 3.01, 
                                 jvar = 3.01)
  expect_true(is.numeric(res1))
  expect_equal(length(res1), 21L)
  expect_identical(names(res1), c("Ygamma", paste0("X", 1:20)))
  # Build a small data.frame by replicate
  dat <- as.data.frame(t(replicate(10, simul_data_UniYX_gamma(totdim = 20, 
                                                              ncomp = 2,
                                                              lvar = 3.01, 
                                                              jvar = 3.01))))
  expect_s3_class(dat, "data.frame")
  expect_equal(ncol(dat), 21L)
  expect_true(all(vapply(dat, is.numeric, logical(1))))
})
