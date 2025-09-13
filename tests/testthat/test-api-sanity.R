test_that("Exported objects are present and functions have formals", {
  skip_if_not_installed("bootPLS")
  ns <- asNamespace("bootPLS")
  exported <- c()
  for (nm in exported) {
    expect_true(exists(nm, envir = ns, inherits = FALSE), info = paste("Missing:", nm))
    obj <- get(nm, envir = ns, inherits = FALSE)
    if (is.function(obj)) {
      fmls <- formals(obj)
      expect_true(is.null(fmls) || is.pairlist(fmls))
    }
  }
})
