####### Finding the number of components using bootstrap for a given eta
#' Title
#'
#' @param x Matrix of predictors.
#' @param y Vector or matrix of responses.
#' @param fold Number of fold for cross-validation
#' @param eta Thresholding parameter. eta should be between 0 and 1.
#' @param R Number of resamplings.
#' @param maxnt Maximum number of components allowed in a spls model.
#' @param kappa Parameter to control the effect of the concavity of the 
#' objective function and the closeness of original and surrogate 
#' direction vectors. kappa is relevant only when responses are multivariate. 
#' kappa should be between 0 and 0.5. Default is 0.5.
#' @param select PLS algorithm for variable selection. Alternatives are 
#' "pls2" or "simpls". Default is "pls2".
#' @param fit PLS algorithm for model fitting. Alternatives are "kernelpls", 
#' "widekernelpls", "simpls", or "oscorespls". Default is "simpls".
#' @param scale.x Scale predictors by dividing each predictor variable by 
#' its sample standard deviation?
#' @param scale.y Scale responses by dividing each response variable by its 
#' sample standard deviation?
#' @param plot.it Plot the results.
#' @param typeBCa Include computation for BCa type interval.
#' @param verbose Displays information on the algorithm.
#'
#' @return list of 3: mspemat matrix of results, eta.opt numeric value, K.opt numeric value)
#' 
#' @export
#'
#' @author Jérémy Magnanensi, Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
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
#' accepted.
#' 
#' @examples
#' set.seed(314)
#' data(pine, package = "plsRglm")
#' Xpine<-pine[,1:10]
#' ypine<-log(pine[,11])
#' nbcomp.bootspls(x=Xpine,y=ypine,eta=.2, maxnt=1)
#' \donttest{
#' set.seed(314)
#' data(pine, package = "plsRglm")
#' Xpine<-pine[,1:10]
#' ypine<-log(pine[,11])
#' nbcomp.bootspls.para(x=Xpine,y=ypine,eta=c(.2,.6))
#' }
nbcomp.bootspls=function (x, y, fold = 10, eta, R=500, maxnt=10, kappa = 0.5, 
                   select = "pls2", fit = "simpls", scale.x = TRUE, 
                   scale.y = FALSE, plot.it = TRUE, typeBCa = TRUE, 
                   verbose=TRUE) 
{
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ip <- c(1:p)
  y <- as.matrix(y)
  q <- ncol(y)
  type <- correctp.withoutK(x, y, eta, kappa, select, fit)
  eta <- type$eta
  kappa <- type$kappa
  select <- type$select
  fit <- type$fit
  foldi <- split(sample(1:n), rep(1:fold, length = n))
  mspemat <- matrix(0, length(eta), 1)
  K.opti=rep(0,length(eta))
  ncolBoot <- 5+2*typeBCa
  for (i in 1:length(eta)) {
    if(verbose){cat(paste("eta =", eta[i], "\n"))}
    ### Finding the number of components
    inter=0
    indK=1
    if(verbose){print(indK)}
    resK=spls.Cboot(x, y, eta = eta[i], kappa = kappa, 
                    K = indK, select = select, fit = fit, scale.x = scale.x, 
                    scale.y = scale.y, verbose = FALSE)
    compTsim=resK$tt
    databoot=cbind(scale(y,scale=F),compTsim)
    bootcomp<-boot::boot(data=databoot, statistic=coefs.plsR.CSim, sim="ordinary", stype="i", R=R)
    confYT=t(as.matrix(plsRglm::confints.bootpls(bootcomp, typeBCa = typeBCa)))
    
    while ((confYT[1,ncolBoot]>0) & (indK<maxnt+1)) {
      indK=indK+1
      if(verbose){print(indK)}
      resK<-spls.Cboot(x, y, eta = eta[i], kappa = kappa, 
                       K = indK, select = select, fit = fit, scale.x = scale.x, 
                       scale.y = scale.y, verbose = FALSE)
      ind=1
      compTsim=resK$tt[,ind]
      databoot=cbind(scale(y,scale=F),compTsim)
      bootcomp<-boot::boot(data=databoot, statistic=coefs.plsR.CSim, sim="ordinary", stype="i", R=R)
      confC=t(as.matrix(plsRglm::confints.bootpls(bootcomp, typeBCa = typeBCa)))
      while ((confC[1,ncolBoot]>0) && (ind<indK)){
        ind=ind+1
        compTsim=resK$tt[,1:ind]
        databoot=cbind(scale(y,scale=F),compTsim)
        bootcomp<-boot::boot(data=databoot, statistic=coefs.plsR.CSim, sim="ordinary", stype="i", R=R)
        confC=t(as.matrix(plsRglm::confints.bootpls(bootcomp, typeBCa = typeBCa)))
      }
      if (ind!=indK){
        #confYT=matrix(0, indK, ncolBoot+1)
        confYT=matrix(0, 1, ncolBoot+1)
      }else{
        confYT=confC
      }
    }
    K.opti[i]=indK - 1
    
    
    mspemati <- matrix(0, fold, 1)
    for (j in 1:fold) {
      omit <- foldi[[j]]
      object <- spls::spls(x[-omit, , drop = FALSE], y[-omit, 
                                                       , drop = FALSE], eta = eta[i], kappa = kappa, 
                           K = K.opti[i], select = select, fit = fit, scale.x = scale.x, 
                           scale.y = scale.y, trace = FALSE)
      newx <- x[omit, , drop = FALSE]
      newx <- scale(newx, object$meanx, object$normx)
      betamat <- object$betamat
      pred <- newx %*% betamat[[K.opti[i]]] + matrix(1, nrow(newx), 
                                                     1) %*% object$mu
      mspemati[j, ] <- mean(apply((y[omit, 
      ] - pred)^2, 2, mean))
    }
    mspemat[i, ] <- apply(mspemati, 2, mean)
  }
  minpmse <- min(mspemat)
  rownames(mspemat) <- paste("eta=",eta,"; K=",K.opti)
  colnames(mspemat) <- ""
  eta.opt <- max(eta[mspemat == minpmse])
  K.opt <- K.opti[eta == max(eta[mspemat == minpmse])]
  if(verbose){cat(paste("\nOptimal parameters: eta = ", eta.opt, ", ", 
            sep = ""))}
  if(verbose){cat(paste("K = ", K.opt, "\n", sep = ""))}
  if (plot.it) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mar=c(5,7,4,2))
    spls::heatmap.spls(t(mspemat), main = "CV MSPE Plot", 
                       coln = 16, as = "n")
  }
  rownames(mspemat) <- paste("eta=", eta,", K=", K.opti)
  cv <- list(mspemat = mspemat, eta.opt = eta.opt, K.opt = K.opt)
  return(cv)
}

####### Finding the number of components using boostrap for a given eta (parallel version)
#' Title
#'
#' @param x Matrix of predictors.
#' @param y Vector or matrix of responses.
#' @param fold Number of fold for cross-validation
#' @param eta Thresholding parameter. eta should be between 0 and 1.
#' @param R Number of resamplings.
#' @param maxnt Maximum number of components allowed in a spls model.
#' @param kappa Parameter to control the effect of the concavity of the 
#' objective function and the closeness of original and surrogate 
#' direction vectors. kappa is relevant only when responses are multivariate. 
#' kappa should be between 0 and 0.5. Default is 0.5.
#' @param select PLS algorithm for variable selection. Alternatives are 
#' "pls2" or "simpls". Default is "pls2".
#' @param fit PLS algorithm for model fitting. Alternatives are "kernelpls", 
#' "widekernelpls", "simpls", or "oscorespls". Default is "simpls".
#' @param scale.x Scale predictors by dividing each predictor variable by 
#' its sample standard deviation?
#' @param scale.y Scale responses by dividing each response variable by its 
#' sample standard deviation?
#' @param plot.it Plot the results.
#' @param typeBCa Include computation for BCa type interval.
#' @param ncpus Number of cpus for parallel computing.
#' @param verbose Displays information on the algorithm.
#'
#' @return list of 3: mspemat matrix of results, eta.opt numeric value, K.opt numeric value)
#' @export
#'
#' @author Jérémy Magnanensi, Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
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
#' accepted.
#' 
#' @examples
#' set.seed(314)
#' data(pine, package = "plsRglm")
#' Xpine<-pine[,1:10]
#' ypine<-log(pine[,11])
#' nbcomp.bootspls.para(x=Xpine,y=ypine,eta=.2, maxnt=1)
#' \donttest{
#' set.seed(314)
#' data(pine, package = "plsRglm")
#' Xpine<-pine[,1:10]
#' ypine<-log(pine[,11])
#' nbcomp.bootspls.para(x=Xpine,y=ypine,eta=c(.2,.6))
#' }
nbcomp.bootspls.para=function (x, y, fold = 10, eta, R=500, maxnt=10, kappa = 0.5, 
                         select = "pls2", fit = "simpls", scale.x = TRUE, 
                         scale.y = FALSE, plot.it = TRUE, typeBCa = TRUE, 
                         ncpus=1, verbose=TRUE) 
{
  lll<-NA
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ip <- c(1:p)
  y <- as.matrix(y)
  q <- ncol(y)
  ncolBoot <- 5+2*typeBCa
  type <- correctp.withoutK(x, y, eta, kappa, select, fit)
  eta <- type$eta
  kappa <- type$kappa
  select <- type$select
  fit <- type$fit
  foldi <- split(sample(1:n), rep(1:fold, length = n))
  #cl <- makeCluster(6)
  #registerDoParallel(cl)
  doParallel::registerDoParallel(cores=ncpus)
  requireNamespace("foreach",quietly = TRUE)
  par.K.opti=foreach::foreach(lll=1:length(eta), .combine="cbind", 
                              .export=c("x","y", "eta", "kappa", "select", 
                                        "fit", "scale.x", "scale.y", "maxnt")) %dopar% {
    if(verbose){print(paste("eta =", eta[lll]))}
    ### Finding the number of components
    inter=0
    indK=1
    #if(verbose){print(indK)}
    resK=spls.Cboot(x, y, eta = eta[lll], kappa = kappa, 
                    K = indK, select = select, fit = fit, scale.x = scale.x, 
                    scale.y = scale.y, verbose = FALSE)
    compTsim=resK$tt
    databoot=cbind(scale(y,scale=F),compTsim)
    bootcomp<-boot::boot(data=databoot, statistic=coefs.plsR.CSim, sim="ordinary", stype="i", R=R)
    confYT=t(as.matrix(plsRglm::confints.bootpls(bootcomp, typeBCa = typeBCa)))
    
    while ((confYT[1,ncolBoot]>0) & (indK<maxnt+1)) {
      indK=indK+1
      if(verbose){print(indK)}
      resK<-spls.Cboot(x, y, eta = eta[lll], kappa = kappa, 
                       K = indK, select = select, fit = fit, scale.x = scale.x, 
                       scale.y = scale.y, verbose = FALSE)
      ind=1
      compTsim=resK$tt[,ind]
      databoot=cbind(scale(y,scale=F),compTsim)
      bootcomp<-boot::boot(data=databoot, statistic=coefs.plsR.CSim, sim="ordinary", stype="i", R=R)
      confC=t(as.matrix(plsRglm::confints.bootpls(bootcomp, typeBCa = typeBCa)))
      while ((confC[1,ncolBoot]>0) && (ind<indK)){
        ind=ind+1
        compTsim=resK$tt[,1:ind]
        databoot=cbind(scale(y,scale=F),compTsim)
        bootcomp<-boot::boot(data=databoot, statistic=coefs.plsR.CSim, sim="ordinary", stype="i", R=R)
        confC=t(as.matrix(plsRglm::confints.bootpls(bootcomp, typeBCa = typeBCa)))
      }
      if (ind!=indK){
        #confYT=matrix(0, indK, ncolBoot+1)
        confYT=matrix(0, 1, ncolBoot+1)
      }else{
        confYT=confC
      }
    }

    mspemati <- matrix(0, fold, 1)
    for (j in 1:fold) {
      omit <- foldi[[j]]
      object <- spls::spls(x[-omit, , drop = FALSE], y[-omit, 
                                                       , drop = FALSE], eta = eta[lll], kappa = kappa, 
                           K = indK - 1, select = select, fit = fit, scale.x = scale.x, 
                           scale.y = scale.y, trace = FALSE)
      newx <- x[omit, , drop = FALSE]
      newx <- scale(newx, object$meanx, object$normx)
      betamat <- object$betamat
      pred <- newx %*% betamat[[indK - 1]] + matrix(1, nrow(newx), 
                                                    1) %*% object$mu
      mspemati[j, ] <- mean(apply((y[omit, 
      ] - pred)^2, 2, mean))
    }
    c(indK - 1, apply(mspemati, 2, mean))
  }
  doParallel::stopImplicitCluster()
  
  if(!is.matrix(par.K.opti)){par.K.opti <- matrix(par.K.opti,ncol=1)}
  K.opti=par.K.opti[1,]
  mspemat=as.matrix(par.K.opti[2,],ncol=1)
  minpmse <- min(mspemat)
  rownames(mspemat) <- paste("eta=",eta,"; K=",K.opti)
  colnames(mspemat) <- ""
  eta.opt <- max(eta[mspemat == minpmse])
  K.opt <- K.opti[eta == max(eta[mspemat == minpmse])]
  if(verbose){cat(paste("\nOptimal parameters: eta = ", eta.opt, ", ", 
            sep = ""))}
  if(verbose){cat(paste("K = ", K.opt, "\n", sep = ""))}
  if (plot.it) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mar=c(5,7,4,2))
    spls::heatmap.spls(t(mspemat), main = "CV MSPE Plot", 
                       coln = 16, as = "n")
  }
  cv <- list(mspemat = mspemat, eta.opt = eta.opt, K.opt = K.opt)
  return(cv)
}