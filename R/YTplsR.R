#' Non-parametric (Y,T) Bootstrap for selecting the number of components in PLSR 
#' models
#'
#' Provides a wrapper for the bootstrap function \code{boot} from the
#' \code{boot} R package.\cr Implements non-parametric bootstraps for PLS
#' Regression models by (Y,T) resampling to select the number of components.
#' 
#' More details on bootstrap techniques are available in the help of the
#' \code{\link[boot:boot]{boot}} function.
#' 

#' @param Y Vector of response.
#' @param X Matrix of predictors.
#' @param R The number of bootstrap replicates. Usually this will be a single
#' positive integer. For importance resampling, some resamples may use one set
#' of weights and others use a different set of weights. In this case \code{R}
#' would be a vector of integers where each component gives the number of
#' resamples from each of the rows of weights.
#' @param sim A character string indicating the type of simulation required.
#' Possible values are \code{"ordinary"} (the default), \code{"balanced"},
#' \code{"permutation"}, or \code{"antithetic"}.
#' @param ncpus	integer: number of processes to be used in parallel operation: 
#' typically one would chose this to the number of available CPUs. 
#' @param parallel The type of parallel operation to be used (if any). 
#' If missing, the default is taken from the option "boot.parallel" 
#' (and if that is not set, "no").
#' @param typeBCa Compute BCa type intervals ?
#' @param verbose Display info during the run of algorithm?
#' 
#' @return A numeric, the number of components selected by the bootstrap.
#'
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
#' data(pine, package="plsRglm")
#' Xpine<-pine[,1:10]
#' ypine<-log(pine[,11])
#' res <- nbcomp.bootplsR(ypine, Xpine)
#' nbcomp.bootplsR(ypine, Xpine, typeBCa=FALSE)
#' \donttest{
#' nbcomp.bootplsR(ypine, Xpine, typeBCa=FALSE, verbose=FALSE)
#' try(nbcomp.bootplsR(ypine, Xpine, sim="permutation"))
#' nbcomp.bootplsR(ypine, Xpine, sim="permutation", typeBCa=FALSE)
#' }
#' 
nbcomp.bootplsR<-function(Y,X,R=500,sim="ordinary",ncpus=1,parallel="no",typeBCa=TRUE, verbose=TRUE){
  indboot2=1
  ncolBoot <- 5+2*typeBCa
  if(verbose){print(indboot2)}
  ressimYT<-plsRglm::PLS_lm(Y, X, nt = indboot2, modele = "pls", scaleX=TRUE, verbose=verbose)
  compTsim=ressimYT$tt
  databoot=cbind(ressimYT$RepY,compTsim)
  if(sim!="permutation"){
    sim.bootSim2<-boot::boot(data=databoot, parallel=parallel, ncpus=ncpus, statistic=coefs.plsR.CSim, sim=sim, stype="i", R=R)
    confYT=t(as.matrix(plsRglm::confints.bootpls(sim.bootSim2,typeBCa = typeBCa)))
    
    while (confYT[1,ncolBoot]>0){
      indboot2=indboot2+1
      if(verbose){print(indboot2)}
      ressimYT<-plsRglm::PLS_lm(Y, X, nt = indboot2, modele = "pls", scaleX=TRUE, verbose=verbose)
      if (ncol(ressimYT$tt)==ressimYT$nt){
        compTsim=ressimYT$tt
        databoot=cbind(ressimYT$RepY,compTsim)
        sim.bootSim2<-boot::boot(data=databoot, parallel=parallel, ncpus=ncpus, statistic=coefs.plsR.CSim, sim=sim, stype="i", R=R)
        confYT=t(as.matrix(plsRglm::confints.bootpls(sim.bootSim2,typeBCa = typeBCa)))
      }
      else{confYT=matrix(0,indboot2,ncolBoot+1)}
    }
  } else {
    sim.bootSim2<-boot::boot(data=databoot, parallel=parallel, ncpus=ncpus, statistic=permcoefs.plsR.CSim, sim="permutation", stype="i", R=R)
    confYT=t(as.matrix(plsRglm::confints.bootpls(sim.bootSim2,typeBCa = typeBCa)))
    
    while (confYT[1,ncolBoot]>0){
      indboot2=indboot2+1
      if(verbose){print(indboot2)}
      ressimYT<-plsRglm::PLS_lm(Y, X, nt = indboot2, modele = "pls", scaleX=TRUE, verbose=verbose)
      if (ncol(ressimYT$tt)==ressimYT$nt){
        compTsim=ressimYT$tt
        databoot=cbind(ressimYT$RepY,compTsim)
        sim.bootSim2<-boot::boot(data=databoot, parallel=parallel, ncpus=ncpus, statistic=permcoefs.plsR.CSim, sim="permutation", stype="i", R=R)
        confYT=t(as.matrix(plsRglm::confints.bootpls(sim.bootSim2,typeBCa = typeBCa)))
      }
      else{confYT=matrix(0,indboot2,ncolBoot+1)}
    }
  }
  if(verbose){print(paste("Optimal number of components: K = ", indboot2-1, sep = ""))}
  return(indboot2-1)
}
