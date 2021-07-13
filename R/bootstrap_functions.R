#' @title Bootstrap (Y,T) functions for PLSR
#'
#' @param dataset Dataset with tt
#' @param i Index for resampling
#'
#' @return Coefficient of the last variable in the linear regression 
#' \code{lm(dataset[i,1] ~ dataset[,-1] - 1)} computed using bootstrap 
#' resampling.
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
#' xran=matrix(rnorm(150),30,5)
#' coefs.plsR.CSim(xran,sample(1:30))
coefs.plsR.CSim <- function(dataset,i)
{
  #x=lm(scale(dataset[i,1],scale=F) ~ dataset[i,-1] - 1)
  x=lm(dataset[i,1] ~ dataset[i,-1] - 1)
  #  x$coefficients[ncol(as.matrix(dataset[,-1]))]
  x$coefficients[ncol(dataset)-1]
}

#' @title Permutation bootstrap (Y,T) function for PLSR
#'
#' @param dataset Dataset with tt
#' @param i Index for resampling
#'
#' @return Coefficient of the last variable in the linear regression 
#' \code{lm(dataset[i,1] ~ dataset[,-1] - 1)} computed using permutation 
#' resampling.
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
#' xran=matrix(rnorm(150),30,5)
#' permcoefs.plsR.CSim(xran,sample(1:30))
permcoefs.plsR.CSim <- function(dataset,i)
{
  #x=lm(scale(dataset[i,1],scale=F) ~ dataset[i,-1] - 1)
  x=lm(dataset[i,1] ~ dataset[,-1] - 1)
#  x$coefficients[ncol(as.matrix(dataset[,-1]))]
  x$coefficients[ncol(dataset)-1]
}


#' @title Bootstrap (Y,T) function for PLSGLR
#'
#' A function passed to \code{boot} to perform bootstrap.
#' 
#' 
#' @param dataRepYtt Dataset with tt components to resample
#' @param ind indices for resampling
#' @param nt number of components to use
#' @param modele type of modele to use, see \link{plsRglm}. Not used, 
#' please specify the family instead.
#' @param family glm family to use, see \link{plsRglm}
#' @param maxcoefvalues maximum values allowed for the estimates of the
#' coefficients to discard those coming from singular bootstrap samples
#' @param ifbootfail value to return if the estimation fails on a bootstrap
#' sample
#' @return estimates on a bootstrap sample or \code{ifbootfail} value if the
#' bootstrap computation fails.
#'
#' @return Numeric vector of the components computed using a bootstrap 
#' resampling.
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
#' library(plsRglm)
#' data(aze_compl, package="plsRglm")
#' Xaze_compl<-aze_compl[,2:34]
#' yaze_compl<-aze_compl$y
#' dataset <- cbind(y=yaze_compl,Xaze_compl)
#' modplsglm <- plsRglm::plsRglm(y~.,data=dataset,4,modele="pls-glm-family",family=binomial)
#' dataRepYtt <- cbind(y = modplsglm$RepY, modplsglm$tt)
#' coefs.plsRglm.CSim(dataRepYtt, sample(1:nrow(dataRepYtt)), 4, 
#' family = binomial, maxcoefvalues=10, ifbootfail=0)
coefs.plsRglm.CSim<-function (dataRepYtt, ind, nt, modele, family = NULL, maxcoefvalues,
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
    return(tempcoefs)
  }
  else {
    return(ifbootfail)
  }
}

#' @title Permutation bootstrap (Y,T) function for PLSGLR
#'
#' A function passed to \code{boot} to perform bootstrap.
#' 
#' 
#' @param dataRepYtt Dataset with tt components to resample
#' @param ind indices for resampling
#' @param nt number of components to use
#' @param modele type of modele to use, see \link{plsRglm}. Not used, 
#' please specify the family instead.
#' @param family glm family to use, see \link{plsRglm}
#' @param maxcoefvalues maximum values allowed for the estimates of the
#' coefficients to discard those coming from singular bootstrap samples
#' @param ifbootfail value to return if the estimation fails on a bootstrap
#' sample
#' @return estimates on a bootstrap sample or \code{ifbootfail} value if the
#' bootstrap computation fails.
#'
#' @return Numeric vector of the components computed using a permutation 
#' resampling.
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
#' library(plsRglm)
#' data(aze_compl, package="plsRglm")
#' Xaze_compl<-aze_compl[,2:34]
#' yaze_compl<-aze_compl$y
#' dataset <- cbind(y=yaze_compl,Xaze_compl)
#' modplsglm <- plsRglm::plsRglm(y~.,data=dataset,4,modele="pls-glm-logistic")
#' dataRepYtt <- cbind(y = modplsglm$RepY, modplsglm$tt)
#' permcoefs.plsRglm.CSim(dataRepYtt, sample(1:nrow(dataRepYtt)), 4, 
#' family = binomial, maxcoefvalues=10, ifbootfail=0)
permcoefs.plsRglm.CSim<-function (dataRepYtt, ind, nt, modele, family = NULL, maxcoefvalues,
                              ifbootfail)
{
  dataRepYb = dataRepYtt[ind, 1]
  Tb = dataRepYtt[, -1]
  tempCb = try(glm(dataRepYb ~ Tb, family = family), silent = TRUE)
  tempcoefs <- tempCb$coefficients[-1]
  Cond <- FALSE
  try(Cond <- is.numeric(tempcoefs) & all(abs(tempcoefs) <
                                            maxcoefvalues), silent = TRUE)
  if (Cond) {                  
    return(tempcoefs)
  }
  else {
    return(ifbootfail)
  }
}
