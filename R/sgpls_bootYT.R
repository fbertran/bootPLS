#' @title Bootstrap (Y,T) function for plsRglm
#'
#' @description A function passed to \code{boot} to perform bootstrap.
#' 
#' 
#' @param dataRepYtt Dataset with tt components to resample
#' @param ind indices for resampling
#' @param nt number of components to use
#' @param modele type of modele to use, see \link[plsRglm]{plsRglm}. Not used, 
#' please specify the family instead.
#' @param family glm family to use, see \link[plsRglm]{plsRglm}
#' @param maxcoefvalues maximum values allowed for the estimates of the
#' coefficients to discard those coming from singular bootstrap samples
#' @param ifbootfail value to return if the estimation fails on a bootstrap
#' sample
#' @return Numeric vector of the components computed using a bootstrap 
#' resampling or \code{ifbootfail} value if the
#' bootstrap computation fails.
#' @export coefs.sgpls.CSim
#'
#' @author Jérémy Magnanensi, Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' 
#' @references A new bootstrap-based stopping criterion in PLS component construction,
#' J. Magnanensi, M. Maumy-Bertrand, N. Meyer and F. Bertrand (2016), in The Multiple Facets of Partial Least Squares and Related Methods, 
#' \doi{10.1007/978-3-319-40643-5_18}\cr
#' 
#' A new universal resample-stable bootstrap-based stopping criterion for PLS component construction,
#' J. Magnanensi, F. Bertrand, M. Maumy-Bertrand and N. Meyer, (2017), Statistics and Computing, 27, 757–774. 
#' \doi{10.1007/s11222-016-9651-4}\cr
#' 
#' New developments in Sparse PLS regression, J. Magnanensi, M. Maumy-Bertrand, 
#' N. Meyer and F. Bertrand, (2021), Frontiers in Applied Mathematics and Statistics, 
#' \doi{10.3389/fams.2021.693126}\cr.
#'
#' @examples
#' set.seed(4619)
#' xran=cbind(rbinom(30,1,.2),matrix(rnorm(150),30,5))
#' coefs.sgpls.CSim(xran, ind=sample(1:nrow(xran)), 
#' maxcoefvalues=1e5, ifbootfail=rep(NA,3))
#' 
coefs.sgpls.CSim<-function (dataRepYtt, ind, nt, modele, family = binomial, maxcoefvalues,
                            ifbootfail)
{
  dataRepYb = dataRepYtt[ind, 1]
  Tb = dataRepYtt[ind, -1]
  
  tempCb = try(glm(dataRepYb ~ Tb, family = family),
               silent = TRUE)
  tempcoefs <- tempCb$coefficients[-1]
  Cond <- FALSE
  try(Cond <- is.numeric(tempcoefs) & all(abs(tempcoefs) <
                                            maxcoefvalues), silent = TRUE)
  if (Cond) {
    return(tempcoefs[nt])
  }
  else {
    return(ifbootfail)
  }
}

#' @title Permutation Bootstrap (Y,T) function for plsRglm
#'
#' @param dataRepYtt Dataset with tt components to resample
#' @param ind indices for resampling
#' @param nt number of components to use
#' @param modele type of modele to use, see \link[plsRglm]{plsRglm}. Not used, 
#' please specify the family instead.
#' @param family glm family to use, see \link[plsRglm]{plsRglm}
#' @param maxcoefvalues maximum values allowed for the estimates of the
#' coefficients to discard those coming from singular bootstrap samples
#' @param ifbootfail value to return if the estimation fails on a bootstrap
#' sample
#' @return Numeric vector of the components computed using a bootstrap 
#' resampling or \code{ifbootfail} value if the
#' bootstrap computation fails.
#' @export
#'
#' @author Jérémy Magnanensi, Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' 
#' @references A new bootstrap-based stopping criterion in PLS component construction,
#' J. Magnanensi, M. Maumy-Bertrand, N. Meyer and F. Bertrand (2016), in The Multiple Facets of Partial Least Squares and Related Methods, 
#' \doi{10.1007/978-3-319-40643-5_18}\cr
#' 
#' A new universal resample-stable bootstrap-based stopping criterion for PLS component construction,
#' J. Magnanensi, F. Bertrand, M. Maumy-Bertrand and N. Meyer, (2017), Statistics and Computing, 27, 757–774. 
#' \doi{10.1007/s11222-016-9651-4}\cr
#' 
#' New developments in Sparse PLS regression, J. Magnanensi, M. Maumy-Bertrand, 
#' N. Meyer and F. Bertrand, (2021), Frontiers in Applied Mathematics and Statistics, 
#' \doi{10.3389/fams.2021.693126}\cr.
#'
#' @examples
#' set.seed(4619)
#' xran=cbind(rbinom(30,1,.2),matrix(rnorm(150),30,5))
#' permcoefs.sgpls.CSim(xran, ind=sample(1:nrow(xran)), maxcoefvalues=1e5, 
#' ifbootfail=rep(NA,3))
#' 
permcoefs.sgpls.CSim<-function (dataRepYtt, ind, nt, modele, family = binomial, maxcoefvalues,
                            ifbootfail)
{
  dataRepYb = dataRepYtt[, 1]
  Tb = dataRepYtt[ind, -1]
  
  tempCb = try(glm(dataRepYb ~ Tb, family = family),
               silent = TRUE)
  tempcoefs <- tempCb$coefficients[-1]
  Cond <- FALSE
  try(Cond <- is.numeric(tempcoefs) & all(abs(tempcoefs) <
                                            maxcoefvalues), silent = TRUE)
  if (Cond) {
    return(tempcoefs[nt])
  }
  else {
    return(ifbootfail)
  }
}


#' @title Number of components for SGPLS using (Y,T) bootstrap
#'
#' @param x Matrix of predictors.
#' @param y Vector or matrix of responses.
#' @param fold Number of fold for cross-validation
#' @param eta Thresholding parameter. eta should be between 0 and 1.
#' @param R Number of resamplings.
#' @param scale.x Scale predictors by dividing each predictor variable by 
#' its sample standard deviation?
#' @param maxnt Maximum number of components allowed in a spls model.
#' @param br Apply Firth's bias reduction procedure?
#' @param ftype Type of Firth's bias reduction procedure. Alternatives are 
#' "iden" (the approximated version) or "hat" (the original version). 
#' Default is "iden".
#' @param plot.it Plot the results.
#' @param typeBCa Include computation for BCa type interval.
#' @param stabvalue A value to hard threshold bootstrap estimates computed from
#' atypical resamplings.
#' @param verbose Additionnal information on the algorithm.
#'
#' @return List of four: error matrix, eta optimal, K optimal and the matrix 
#' of results.
#' @export
#'
#' @author Jérémy Magnanensi, Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' 
#' @references A new bootstrap-based stopping criterion in PLS component construction,
#' J. Magnanensi, M. Maumy-Bertrand, N. Meyer and F. Bertrand (2016), in The Multiple Facets of Partial Least Squares and Related Methods, 
#' \doi{10.1007/978-3-319-40643-5_18}\cr
#' 
#' A new universal resample-stable bootstrap-based stopping criterion for PLS component construction,
#' J. Magnanensi, F. Bertrand, M. Maumy-Bertrand and N. Meyer, (2017), Statistics and Computing, 27, 757–774. 
#' \doi{10.1007/s11222-016-9651-4}\cr
#' 
#' New developments in Sparse PLS regression, J. Magnanensi, M. Maumy-Bertrand, 
#' N. Meyer and F. Bertrand, (2021), Frontiers in Applied Mathematics and Statistics, 
#' \doi{10.3389/fams.2021.693126}\cr.
#' 
#' @examples
#' set.seed(4619)
#' data(prostate, package="spls")
#' nbcomp.bootsgpls((prostate$x)[,1:30], prostate$y, R=250, eta=0.2, maxnt=1, typeBCa = FALSE)
#' \donttest{
#' set.seed(4619)
#' data(prostate, package="spls")
#' nbcomp.bootsgpls(prostate$x, prostate$y, R=250, eta=c(0.2,0.6), typeBCa = FALSE)
#' }
nbcomp.bootsgpls=function (x, y, fold = 10, eta, R, scale.x = TRUE, maxnt=10, plot.it = TRUE, 
                           br = TRUE, ftype = "iden", typeBCa=TRUE, stabvalue = 1e+06, verbose=TRUE) 
{
  ifbootfail <- as.matrix(as.numeric(ifelse(any(class(y) == "factor"), rep(NA, ncol(x) + nlevels(y)), rep(NA, ncol(x)+1))))
  result.mat <- c()
  foldi <- cv.split(y, fold)
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ip <- c(1:p)
  y <- as.matrix(y)
  q <- ncol(y)
  K.opti=rep(0,length(eta))
  ncolBoot <- 5+2*typeBCa
  for (lll in 1:length(eta)) {
    if(verbose){print(paste("eta =", eta[lll]))}
    
    ### Finding the number of components
    inter=0
    indK=1
    resK=sgpls.T(x, y, K = indK, eta = eta[lll], scale.x = scale.x, eps = 1e-05, denom.eps = 1e-20, 
                 zero.eps = 1e-05, maxstep = 100, br = br, ftype = ftype)
    if(verbose){print(paste("K =", indK))}
    maxcoefvalues <- stabvalue * abs(resK$CoeffC)
    compTsim=resK$tt
    databoot=cbind(y,compTsim)
    bootcomp<-boot::boot(databoot, statistic=coefs.sgpls.CSim, nt=indK, maxcoefvalues = maxcoefvalues, ifbootfail = ifbootfail, sim="ordinary", stype="i", R=R)
    confYT=t(as.matrix(plsRglm::confints.bootpls(bootcomp, typeBCa = typeBCa)))
    
    while ((confYT[1,ncolBoot]>0) & (indK<maxnt+1)) {
      indK=indK+1             
      resK<-sgpls.T(x, y, K = indK, eta = eta[lll], scale.x = scale.x, eps = 1e-05, denom.eps = 1e-20, 
                    zero.eps = 1e-05, maxstep = 100, br = br, ftype = ftype)
      if(verbose){print(paste("K =", indK))}
      maxcoefvalues <- stabvalue * abs(resK$CoeffC)
      ind=1
      compTsim=resK$tt[,ind]
      databoot=cbind(y,compTsim)
      bootcomp<-boot::boot(databoot, statistic=coefs.sgpls.CSim, nt=ind, maxcoefvalues = maxcoefvalues[ind], ifbootfail = ifbootfail, sim="ordinary", stype="i", R=R)
      confC=t(as.matrix(plsRglm::confints.bootpls(bootcomp, typeBCa = typeBCa)))
      while ((confC[1,ncolBoot]>0) & (ind<indK)){
        ind=ind+1
        compTsim=resK$tt[,1:ind]
        databoot=cbind(y,compTsim)
        bootcomp<-boot::boot(databoot, statistic=coefs.sgpls.CSim, nt=ind, maxcoefvalues = maxcoefvalues[1:ind], ifbootfail = ifbootfail, sim="ordinary", stype="i", R=R)
        confC=t(as.matrix(plsRglm::confints.bootpls(bootcomp, typeBCa = typeBCa)))
      }
      if (ind!=indK){
        confYT=matrix(0, 1, ncolBoot+1)
      }else{
        confYT=confC
      }
    }
    K.opti[lll]=indK - 1
  }
  
  .fit.sgpls.binary <- function(eta.val, K.val) {
    mspemati <- rep(0, fold)
    Ai <- rep(0, fold)
    for (k in 1:fold) {
      omit <- foldi[[k]]
      train.x <- x[-omit, ]
      train.y <- y[-omit, ]
      test.x <- x[omit, ]
      test.y <- y[omit, ]
      sgpls.fit <- spls::sgpls(train.x, train.y, K = K.val, eta = eta.val, 
                               scale.x = scale.x, br = br, ftype = ftype)
      pred <- as.numeric(as.vector(predict(sgpls.fit, newx = test.x)))
      mspemati[k] <- mean(as.numeric(pred != test.y))
      Ai[k] <- mean(length(sgpls.fit$A))
    }
    mspe.ij <- c(mean(mspemati), mean(Ai), eta.val, K.val)
    return(mspe.ij)
  }
  result.mat=matrix(0,length(eta),4)
  for (i in 1:length(eta)){
    result.mat[i,] <- .fit.sgpls.binary(eta.val = eta[i],K.val = K.opti[i]) 
  }                                                                
  mspemat <- matrix(result.mat[, 1], length(eta),1)
  rownames(mspemat) <- paste("eta=",eta,"; K=",K.opti)
  colnames(mspemat) <- ""
  cands <- result.mat[result.mat[, 1] == min(result.mat[, 1]), 
                      , drop = FALSE]
  cands <- cands[cands[, 2] == min(cands[, 2]), , drop = FALSE]
  cands <- cands[cands[, 4] == min(cands[, 4]), , drop = FALSE]
  cands <- cands[cands[, 3] == max(cands[, 3]), , drop = FALSE]
  K.opt <- cands[, 4]
  eta.opt <- cands[, 3]
  if(verbose){cat(paste("\nOptimal parameters: eta = ", eta.opt, ", ", 
                        sep = ""))}
  if(verbose){cat(paste("K = ", K.opt, "\n", sep = ""))}
  if (plot.it) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mar=c(5,7,4,2))
    spls::heatmap.spls(t(mspemat), main = "CV error Plot", 
                       coln = 16, as = "n")
  }
  rownames(mspemat) <- paste("eta=",eta,"; K=",K.opti)
  colnames(mspemat) <- ""
  cv <- list(err.mat = mspemat, eta.opt = eta.opt, K.opt = K.opt, cands=cands)
  return(cv)
}


#' @title Number of components for SGPLS using (Y,T) bootstrap 
#' (parallel version)
#'
#' @param x Matrix of predictors.
#' @param y Vector or matrix of responses.
#' @param fold Number of fold for cross-validation.
#' @param eta Thresholding parameter. eta should be between 0 and 1.
#' @param R Number of resamplings.
#' @param scale.x Scale predictors by dividing each predictor variable by 
#' its sample standard deviation?
#' @param maxnt Maximum number of components allowed in a spls model.
#' @param br Apply Firth's bias reduction procedure?
#' @param ftype Type of Firth's bias reduction procedure. Alternatives are 
#' "iden" (the approximated version) or "hat" (the original version). 
#' Default is "iden".
#' @param ncpus Number of cpus for parallel computing.
#' @param plot.it Plot the results.
#' @param typeBCa Include computation for BCa type interval.
#' @param stabvalue A value to hard threshold bootstrap estimates computed from
#' atypical resamplings.
#' @param verbose Additionnal information on the algorithm.
#'
#' @return List of four: error matrix, eta optimal, K optimal and the matrix 
#' of results.
#' @export
#'
#' @author Jérémy Magnanensi, Frédéric Bertrand\cr
#' \email{frederic.bertrand@@lecnam.net}\cr
#' \url{https://fbertran.github.io/homepage/}
#' 
#' @references A new bootstrap-based stopping criterion in PLS component construction,
#' J. Magnanensi, M. Maumy-Bertrand, N. Meyer and F. Bertrand (2016), in The Multiple Facets of Partial Least Squares and Related Methods, 
#' \doi{10.1007/978-3-319-40643-5_18}\cr
#' 
#' A new universal resample-stable bootstrap-based stopping criterion for PLS component construction,
#' J. Magnanensi, F. Bertrand, M. Maumy-Bertrand and N. Meyer, (2017), Statistics and Computing, 27, 757–774. 
#' \doi{10.1007/s11222-016-9651-4}\cr
#' 
#' New developments in Sparse PLS regression, J. Magnanensi, M. Maumy-Bertrand, 
#' N. Meyer and F. Bertrand, (2021), Frontiers in Applied Mathematics and Statistics, 
#' \doi{10.3389/fams.2021.693126}\cr.
#' 
#' @examples
#' set.seed(4619)
#' data(prostate, package="spls")
#' nbcomp.bootsgpls.para((prostate$x)[,1:30], prostate$y, R=250, eta=0.2, maxnt=1, typeBCa = FALSE)
#' \donttest{
#' set.seed(4619)
#' data(prostate, package="spls")
#' nbcomp.bootsgpls.para(prostate$x, prostate$y, R=250, eta=c(0.2,0.6), typeBCa = FALSE)
#' }
nbcomp.bootsgpls.para=function (x, y, fold = 10, eta, R, scale.x = TRUE, maxnt=10, 
                                br = TRUE, ftype = "iden", ncpus = 1, plot.it = TRUE, 
                                typeBCa=TRUE, stabvalue = 1e+06, verbose=TRUE) 
{
  lll<-NA
  ifbootfail <- as.matrix(as.numeric(ifelse(any(class(y) == "factor"), rep(NA, ncol(x) + nlevels(y)), rep(NA, ncol(x)+1))))
  result.mat <- c()
  foldi <- cv.split(y, fold)
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ip <- c(1:p)
  y <- as.matrix(y)
  q <- ncol(y)
  K.opti=rep(0,length(eta))
  ncolBoot <- 5+2*typeBCa
  
  doParallel::registerDoParallel(cores=ncpus)
  requireNamespace("foreach",quietly = TRUE)
  par.K.opti=foreach::foreach(lll=1:length(eta), .combine="cbind", 
                              .export=c("eta", "ftype", "br", "scale.x", 
                                        "R", "maxnt", "stabvalue")) %dopar% {
                                          if(verbose){print(paste("eta =", eta[lll]))}
                                          
                                          ### Finding the number of components
                                          inter=0
                                          indK=1
                                          resK=sgpls.T(x, y, K = indK, eta = eta[lll], scale.x = scale.x, eps = 1e-05, denom.eps = 1e-20, 
                                                       zero.eps = 1e-05, maxstep = 100, br = br, ftype = ftype)
                                          if(verbose){print(paste("K =", indK))}
                                          maxcoefvalues <- stabvalue * abs(resK$CoeffC)
                                          compTsim=resK$tt
                                          databoot=cbind(y,compTsim)
                                          bootcomp<-boot::boot(databoot, statistic=coefs.sgpls.CSim, nt=indK, maxcoefvalues = maxcoefvalues, ifbootfail = ifbootfail, sim="ordinary", stype="i", R=R)
                                          confYT=t(as.matrix(plsRglm::confints.bootpls(bootcomp, typeBCa = typeBCa)))
                                          
                                          while ((confYT[1,ncolBoot]>0) & (indK<maxnt+1)) {
                                            indK=indK+1             
                                            resK<-sgpls.T(x, y, K = indK, eta = eta[lll], scale.x = scale.x, eps = 1e-05, denom.eps = 1e-20, 
                                                          zero.eps = 1e-05, maxstep = 100, br = br, ftype = ftype)
                                            if(verbose){print(paste("K =", indK))}
                                            maxcoefvalues <- stabvalue * abs(resK$CoeffC)
                                            ind=1
                                            compTsim=resK$tt[,ind]
                                            databoot=cbind(y,compTsim)
                                            bootcomp<-boot::boot(databoot, statistic=coefs.sgpls.CSim, nt=ind, maxcoefvalues = maxcoefvalues[ind], ifbootfail = ifbootfail, sim="ordinary", stype="i", R=R)
                                            confC=t(as.matrix(plsRglm::confints.bootpls(bootcomp, typeBCa = typeBCa)))
                                            while ((confC[1,ncolBoot]>0) & (ind<indK)){
                                              ind=ind+1
                                              compTsim=resK$tt[,1:ind]
                                              databoot=cbind(y,compTsim)
                                              bootcomp<-boot::boot(databoot, statistic=coefs.sgpls.CSim, nt=ind, maxcoefvalues = maxcoefvalues[1:ind], ifbootfail = ifbootfail, sim="ordinary", stype="i", R=R)
                                              confC=t(as.matrix(plsRglm::confints.bootpls(bootcomp, typeBCa = typeBCa)))
                                            }
                                            if (ind!=indK){
                                              confYT=matrix(0, 1, ncolBoot+1)
                                            }else{
                                              confYT=confC
                                            }
                                          }
                                          K.opti[lll]=indK - 1
                                        }
  doParallel::stopImplicitCluster()
  
  .fit.sgpls.binary <- function(eta.val, K.val) {
    mspemati <- rep(0, fold)
    Ai <- rep(0, fold)
    for (k in 1:fold) {
      omit <- foldi[[k]]
      train.x <- x[-omit, ]
      train.y <- y[-omit, ]
      test.x <- x[omit, ]
      test.y <- y[omit, ]
      sgpls.fit <- spls::sgpls(train.x, train.y, K = K.val, eta = eta.val, 
                               scale.x = scale.x, br = br, ftype = ftype)
      pred <- as.numeric(as.vector(predict(sgpls.fit, newx = test.x)))
      mspemati[k] <- mean(as.numeric(pred != test.y))
      Ai[k] <- mean(length(sgpls.fit$A))
    }
    mspe.ij <- c(mean(mspemati), mean(Ai), eta.val, K.val)
    return(mspe.ij)
  }
  result.mat=matrix(0,length(eta),4)
  for (i in 1:length(eta)){
    result.mat[i,] <- .fit.sgpls.binary(eta.val = eta[i],K.val = par.K.opti[i]) 
  }                                                                
  mspemat <- matrix(result.mat[, 1], length(eta),1)
  rownames(mspemat) <- paste("eta=",eta,"; K=",par.K.opti)
  colnames(mspemat) <- ""
  cands <- result.mat[result.mat[, 1] == min(result.mat[, 1]), 
                      , drop = FALSE]
  cands <- cands[cands[, 2] == min(cands[, 2]), , drop = FALSE]
  cands <- cands[cands[, 4] == min(cands[, 4]), , drop = FALSE]
  cands <- cands[cands[, 3] == max(cands[, 3]), , drop = FALSE]
  K.opt <- cands[, 4]
  eta.opt <- cands[, 3]
  if(verbose){cat(paste("\nOptimal parameters: eta = ", eta.opt, ", ", 
                        sep = ""))}
  if(verbose){cat(paste("K = ", K.opt, "\n", sep = ""))}
  if (plot.it) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mar=c(5,7,4,2))
    spls::heatmap.spls(t(mspemat), main = "CV error Plot", 
                       coln = 16, as = "n")
  }
  rownames(mspemat) <- paste("eta=",eta,"; K=",par.K.opti)
  colnames(mspemat) <- ""
  cv <- list(err.mat = mspemat, eta.opt = eta.opt, K.opt = K.opt, cands=cands)
  return(cv)
}
