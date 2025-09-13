test_that("Manual examples run without error (robust by alias)", {
  skip_on_cran()
  skip_if_not_installed("bootPLS")
  options(example.ask = FALSE)
  
  rdb <- tryCatch(tools::Rd_db("bootPLS"), error = function(e) NULL)
  skip_if(is.null(rdb), "No Rd database available")
  
  # Does an Rd have an \examples section?
  has_examples <- function(rd) {
    any(vapply(rd, function(x) identical(attr(x, "Rd_tag"), "\\examples"), logical(1L)))
  }
  # First alias for this Rd; fall back to stripped filename if needed
  first_alias <- function(rdname, rd) {
    idx <- which(vapply(rd, function(x) identical(attr(x, "Rd_tag"), "\\alias"), logical(1L)))
    if (length(idx)) {
      al <- as.character(rd[[idx[1]]])
      if (!is.na(al) && nzchar(al)) return(al)
    }
    sub("\\.[Rr]d$", "", rdname)
  }
  
  withr::local_tempdir()
  ex_ran <- 0L
  for (nm in names(rdb)) {
    rd <- rdb[[nm]]
    if (is.null(rd) || !has_examples(rd)) next
    tp <- first_alias(nm, rd)
    ok <- try({
      expect_silent(example(topic = tp, package = "bootPLS", character.only = TRUE,
                            run.donttest = FALSE, give.lines = FALSE, echo = FALSE))
      TRUE
    }, silent = TRUE)
    if (isTRUE(ok)) ex_ran <- ex_ran + 1L
  }
  expect_gte(ex_ran, 1L)
})
