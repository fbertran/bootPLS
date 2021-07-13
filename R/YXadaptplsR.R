#' @title Bootstrap (Y,X) for the coefficients with number of components updated 
#' for each resampling.
#'
#' @param dataset Dataset to use.
#' @param i Vector of resampling.
#' @param R Number of resamplings to find the number of components.
#' @param ncpus	integer: number of processes to be used in parallel operation: 
#' typically one would chose this to the number of available CPUs.
#' @param parallel The type of parallel operation to be used (if any). 
#' If missing, the default is taken from the option "boot.parallel" 
#' (and if that is not set, "no").
#' @param verbose Suppress information messages.
#'
#' @return Numeric vector: first value is the number of components, the 
#' remaining values are the coefficients the variables computed for that number 
#' of components. 
#' @export
#'
#' @author Jérémy Magnanensi, Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{http://github.fbertran.io/homepage/}
#' 
#' @references A new bootstrap-based stopping criterion in PLS component construction,
#' J. Magnanensi, M. Maumy-Bertrand, N. Meyer and F. Bertrand (2016), in The Multiple Facets of Partial Least Squares and Related Methods, 
#' \doi{10.1007/978-3-319-40643-5_18}\cr
#' 
#' A new universal resample-stable bootstrap-based stopping criterion for PLS component construction,
#' J. Magnanensi, F. Bertrand, M. Maumy-Bertrand and N. Meyer, (2017), Statistics and Compututing, 27, 757–774. 
#' \doi{10.1007/s11222-016-9651-4}\cr
#' 
#' New developments in Sparse PLS regression, J. Magnanensi, M. Maumy-Bertrand, 
#' N. Meyer and F. Bertrand, (2021), Frontiers in Applied Mathematics and Statistics, 
#' \doi{10.3389/fams.2021.693126}
#' 
#' @examples
#' set.seed(314)
#' ncol=10
#' xran=matrix(rnorm(30*ncol),30,ncol)
#' coefs.plsR.adapt.ncomp(xran,sample(1:30))
#' \donttest{
#' coefs.plsR.adapt.ncomp(xran,sample(1:30),ncpus=2,parallel="multicore")
#' }
coefs.plsR.adapt.ncomp <- function(dataset,i,R=1000,ncpus=1,parallel="no",verbose=FALSE)          
{
  nbcomp=bootPLS::nbcomp.bootplsR(dataset[i,1],dataset[i,-1],R=R,ncpus=ncpus,parallel=parallel, verbose =FALSE)
  if (nbcomp!=0){
    ressim<-plsRglm::PLS_glm(dataset[i,1],dataset[i,-1], nt = nbcomp, modele = "pls", scaleX=TRUE, verbose =FALSE)
    c(nbcomp,ressim$Coeffs)
  }
  else {c(0,mean(dataset[i,1]),rep(0,ncol(dataset)-1))}
}


