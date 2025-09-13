test_that("Exported objects exist in namespace", {
  exported <- c('coefs.plsR.CSim', 'coefs.plsR.adapt.ncomp', 'coefs.plsRglm.CSim', 'coefs.sgpls.CSim', 'nbcomp.bootplsR', 'nbcomp.bootplsRglm', 'nbcomp.bootsgpls', 'nbcomp.bootsgpls.para', 'nbcomp.bootspls', 'nbcomp.bootspls.para', 'permcoefs.plsR.CSim', 'permcoefs.plsRglm.CSim', 'permcoefs.sgpls.CSim', 'signpred2', 'simul_data_UniYX_gamma')
  for (nm in exported) {
    expect_true(exists(nm, envir = asNamespace("bootPLS"), inherits = FALSE), info = paste("Missing:", nm))
  }
})
