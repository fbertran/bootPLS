test_that("All example datasets load and have sane structure", {
  skip_if_not_installed("bootPLS")
  ds <- data(package = "bootPLS")$results
  if (is.null(ds) || nrow(ds) == 0L) skip("No packaged datasets found")
  n_ok <- 0L
  for (nm in ds[, "Item"]) {
    env <- new.env(parent = emptyenv())
    expect_silent(data(list = nm, package = "bootPLS", envir = env))
    expect_true(exists(nm, envir = env), info = paste("Dataset missing:", nm))
    obj <- get(nm, envir = env)
    # Basic structural expectations
    expect_true(is.data.frame(obj) || is.matrix(obj) || is.numeric(obj) || is.list(obj),
                info = paste("Unexpected dataset type for", nm))
    n_ok <- n_ok + 1L
  }
  expect_gte(n_ok, 1L)
})
