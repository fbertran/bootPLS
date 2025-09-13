test_that("inst/extdata files (if any) are readable", {
  skip_if_not_installed("bootPLS")
  ext <- system.file("extdata", package = "bootPLS")
  if (identical(ext, "") || !dir.exists(ext)) skip("No inst/extdata")
  files <- list.files(ext, full.names = TRUE)
  if (!length(files)) skip("inst/extdata present but empty")
  for (fp in files) {
    if (grepl("\\.(csv|tsv|txt)$", fp, ignore.case = TRUE)) {
      sep <- if (grepl("\\.tsv$", fp, ignore.case = TRUE)) "\t" else ","
      dat <- tryCatch(utils::read.table(fp, header = TRUE, sep = sep, quote = '"', comment.char = ""),
                      error = function(e) NULL)
      expect_false(is.null(dat), info = paste("Failed to read", basename(fp)))
      expect_gt(NROW(dat), 0, info = paste("Empty data in", basename(fp)))
    } else {
      # Smoke check: file must exist and be non-zero
      expect_true(file.exists(fp))
      info <- file.info(fp)
      expect_gt(info$size, 0, info = paste("Zero-sized file:", basename(fp)))
    }
  }
})
