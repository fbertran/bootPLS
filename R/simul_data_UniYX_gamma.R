#' Data generating function for univariate gamma plsR models
#' 
#' This function generates a single univariate gamma response value \eqn{Ygamma}
#' and a vector of explanatory variables \eqn{(X_1,\ldots,X_{totdim})} drawn
#' from a model with a given number of latent components.
#' 
#' This function should be combined with the replicate function to give rise to
#' a larger dataset. The algorithm used is a modification of a port of the one
#' described in the article of Li which is a multivariate generalization of the
#' algorithm of Naes and Martens.
#' 
#' @param totdim Number of columns of the X vector (from \code{ncomp} to
#' hardware limits)
#' @param ncomp Number of latent components in the model (to use noise, select ncomp=3)
#' @param jvar First variance parameter
#' @param lvar Second variance parameter
#' @param link Character specification of the link function in the mean model
#' (mu). Currently, "\code{inverse}", "\code{log}" and "\code{identity}" are supported.
#' Alternatively, an object of class "link-glm" can be supplied.
#' @param offset Offset on the linear scale
#' @return \item{vector}{\eqn{(Ygamma,X_1,\ldots,X_{totdim})}}
#' @author Jeremy Magnanensi, Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso \code{\link[plsRglm]{simul_data_UniYX}}
#' @references T. Naes, H. Martens, Comparison of prediction methods for
#' multicollinear data, Commun. Stat., Simul. 14 (1985) 545-576.\cr
#' %\url{http://dx.doi.org/10.1080/03610918508812458}\cr Baibing Li, Julian
#' Morris, Elaine B. Martin, Model selection for partial least squares
#' regression, Chemometrics and Intelligent Laboratory Systems 64 (2002), 
#' 79-89, \doi{10.1016/S0169-7439(02)00051-5}.
#' @keywords datagen utilities
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
#' \doi{10.3389/fams.2021.693126}\cr.
#' 
#' @export simul_data_UniYX_gamma
#' @examples
#' set.seed(314)
#' ncomp=rep(3,100)
#' totdimpos=7:50
#' totdim=sample(totdimpos,100,replace=TRUE)
#' l=3.01
#' #for (l in seq(3.01,15.51,by=0.5)) {
#' j=3.01
#' #for (j in seq(3.01,9.51,by=0.5))  {
#' i=44
#' #for ( i in 1:100){
#' set.seed(i)
#' totdimi<-totdim[i]
#' ncompi<-ncomp[i]
#' datasim <- t(replicate(200,simul_data_UniYX_gamma(totdimi,ncompi,j,l)))
#' #}
#' #}
#' #}
#' pairs(datasim)
#' rm(i,j,l,totdimi,ncompi,datasim)
simul_data_UniYX_gamma=function (totdim, ncomp, jvar, lvar, link = "inverse", offset = 0)
{
  if (is.character(link)) {
    linkstr <- link
    linkobj <- make.link(linkstr)
  }
  else {
    linkobj <- link
    linkstr <- link$name
  }
  dimok <- FALSE
  varsR <- c(6, 5.5, 5, jvar, lvar, 1/2)
  varepsilon <- 0.01
  varsF <- c(12.25, 12.125, 12.05, 12.0125, 12.005, 12.00125)
  if (totdim == 1) {
    stop("'totdim' must be > 1")
  }
  if (ncomp == 1) {
    stop("'ncomp' must be > 1")
  }
  if (totdim == 2) {
    dimok <- TRUE
    if (!(ncomp %in% 1:totdim)) {
      warning(paste("ncomp must be <= ", totdim, "\n"))
      warning(paste("ncomp was set to ", totdim, "\n"))
      ncomp <- totdim
    }
    ksi1 <- c(1, 1)/sqrt(2)
    ksi2 <- c(1, -1)/sqrt(2)
    ksi <- cbind(ksi1, ksi2)[1:totdim, 1:ncomp]
  }
  if (totdim == 3) {
    dimok <- TRUE
    if (!(ncomp %in% 1:totdim)) {
      warning(paste("ncomp must be <= ", totdim, "\n"))
      warning(paste("ncomp was set to ", totdim, "\n"))
      ncomp <- totdim
    }
    ksi1 <- c(1, 1, 1)/sqrt(3)
    ksi2 <- c(-1/2, -1/2, 1)/sqrt(3/2)
    ksi3 <- c(-1, 1, 0)/sqrt(2)
    ksi <- cbind(ksi1, ksi2, ksi3)[1:totdim, 1:ncomp]
  }
  if (totdim == 4) {
    dimok <- TRUE
    if (!(ncomp %in% 1:totdim)) {
      warning(paste("ncomp must be <= ", totdim, "\n"))
      warning(paste("ncomp was set to ", totdim, "\n"))
      ncomp <- totdim
    }
    ksi1 <- c(1, 1, 1, 1)/2
    ksi2 <- c(1, -1, 1, -1)/2
    ksi3 <- c(1, 1, -1, -1)/2
    ksi4 <- c(1, -1, -1, 1)/2
    ksi <- cbind(ksi1, ksi2, ksi3, ksi4)[1:totdim, 1:ncomp]
  }
  if (totdim == 5) {
    dimok <- TRUE
    if (!(ncomp %in% 1:totdim)) {
      warning(paste("ncomp must be <= ", totdim, "\n"))
      warning(paste("ncomp was set to ", totdim, "\n"))
      ncomp <- totdim
    }
    ksi1 <- c(1, 1, 1, 1, 1)/sqrt(5)
    ksi2 <- c(1, -1, 1, -1, 0)/2
    ksi3 <- c(1, 1, -1, -1, 0)/2
    ksi4 <- c(1, -1, -1, 1, 0)/2
    ksi5 <- c(1/4, 1/4, 1/4, 1/4, -1)/sqrt(5) * 2
    ksi <- cbind(ksi1, ksi2, ksi3, ksi4, ksi5)[1:totdim,
                                               1:ncomp]
  }
  if ((totdim%%6 == 0) & (totdim >= 6)) {
    dimok <- TRUE
    if (!(ncomp %in% 1:6)) {
      warning(paste("ncomp must be <= ", 6, "\n"))
      warning(paste("ncomp was set to ", 6, "\n"))
      ncomp <- 6
    }
    ksi1 <- rep(c(1, 1, 1, 1, 1, 1), totdim/6)/sqrt(totdim)
    ksi2 <- rep(c(-1/2, -1/2, 1, -1/2, -1/2, 1), totdim/6)/sqrt(totdim/2)
    ksi3 <- rep(c(-1, 1, 0, -1, 1, 0), totdim/6)/sqrt(2 *
                                                        totdim/3)
    ksi4 <- rep(c(1, 1, 1, -1, -1, -1), totdim/6)/sqrt(totdim)
    ksi5 <- rep(c(-1/2, -1/2, 1, 1/2, 1/2, -1), totdim/6)/sqrt(totdim/2)
    ksi6 <- rep(c(-1, 1, 0, 1, -1, 0), totdim/6)/sqrt(2 *
                                                        totdim/3)
    ksi <- cbind(ksi1, ksi2, ksi3, ksi4, ksi5, ksi6)[1:totdim,
                                                     1:(ncomp+1)]
  }
  if ((totdim%%6 == 1) & (totdim >= 6)) {
    dimok <- TRUE
    if (!(ncomp %in% 1:6)) {
      warning(paste("ncomp must be <= ", 6, "\n"))
      warning(paste("ncomp was set to ", 6, "\n"))
      ncomp <- 6
    }
    ksi1 <- c(1, 1, 1, 1, 1, 1, 1, rep(c(1, 1, 1, 1, 1, 1),
                                       (totdim - 1)/6 - 1))/sqrt(totdim)
    ksi2 <- c(-1/2, -1/2, 1, 1, -1, 1, -1, rep(c(-1/2, -1/2,
                                                 1, -1/2, -1/2, 1), (totdim - 1)/6 - 1))/sqrt(3 *
                                                                                                ((totdim - 1)/6 - 1) + 11/2)
    ksi3 <- c(-1, 1, 0, 1, 1, -1, -1, rep(c(-1, 1, 0, -1,
                                            1, 0), (totdim - 1)/6 - 1))/sqrt(4 * ((totdim - 1)/6 -
                                                                                    1) + 6)
    ksi4 <- c(-8/11, -8/11, 16/11, -6/11, 6/11, -6/11, 6/11,
              rep(c(1, 1, 1, -1, -1, -1), (totdim - 1)/6 - 1))/sqrt(6 *
                                                                      ((totdim - 1)/6 - 1) + 4 * sqrt(33)/11)
    ksi5 <- c(2/3, -2/3, 0, 1/3, 1/3, -1/3, -1/3, rep(c(-1/2,
                                                        -1/2, 1, 1/2, 1/2, -1), (totdim - 1)/6 - 1))/sqrt(3 *
                                                                                                            ((totdim - 1)/6 - 1) + 2 * sqrt(3)/3)
    ksi6 <- c(0, 0, 0, 1, -1, -1, 1, rep(c(-1, 1, 0, 1, -1,
                                           0), (totdim - 1)/6 - 1))/sqrt(4 * ((totdim - 1)/6 -
                                                                                1) + 4)
    ksi <- cbind(ksi1, ksi2, ksi3, ksi4, ksi5, ksi6)[1:totdim,
                                                     1:(ncomp+1)]
  }
  if ((totdim%%6 == 2) & (totdim >= 6)) {
    dimok <- TRUE
    if (!(ncomp %in% 1:6)) {
      warning(paste("ncomp must be <= ", 6, "\n"))
      warning(paste("ncomp was set to ", 6, "\n"))
      ncomp <- 6
    }
    ksi1 <- c(1, 1, 1, 1, 1, 1, 1, 1, rep(c(1, 1, 1, 1, 1,
                                            1), (totdim - 2)/6 - 1))/sqrt(totdim)
    ksi2 <- c(1, -1, 1, -1, 1, -1, 1, -1, rep(c(-1/2, -1/2,
                                                1, -1/2, -1/2, 1), (totdim - 2)/6 - 1))/sqrt(3 *
                                                                                               ((totdim - 2)/6 - 1) + 8)
    ksi3 <- c(1, 1, -1, -1, 1, 1, -1, -1, rep(c(-1, 1, 0,
                                                -1, 1, 0), (totdim - 2)/6 - 1))/sqrt(4 * ((totdim -
                                                                                             2)/6 - 1) + 8)
    ksi4 <- c(1, -1, -1, 1, 1, -1, -1, 1, rep(c(1, 1, 1,
                                                -1, -1, -1), (totdim - 2)/6 - 1))/sqrt(6 * ((totdim -
                                                                                               2)/6 - 1) + 8)
    ksi5 <- c(1, 1, 1, 1, -1, -1, -1, -1, rep(c(-1/2, -1/2,
                                                1, 1/2, 1/2, -1), (totdim - 2)/6 - 1))/sqrt(3 * ((totdim -
                                                                                                    2)/6 - 1) + 8)
    ksi6 <- c(1, -1, 1, -1, -1, 1, -1, 1, rep(c(-1, 1, 0,
                                                1, -1, 0), (totdim - 2)/6 - 1))/sqrt(4 * ((totdim -
                                                                                             2)/6 - 1) + 8)
    ksi <- cbind(ksi1, ksi2, ksi3, ksi4, ksi5, ksi6)[1:totdim,
                                                     1:(ncomp+1)]
  }
  if ((totdim%%6 == 3) & (totdim >= 6)) {
    dimok <- TRUE
    if (!(ncomp %in% 1:6)) {
      warning(paste("ncomp must be <= ", 6, "\n"))
      warning(paste("ncomp was set to ", 6, "\n"))
      ncomp <- 6
    }
    ksi1 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, rep(c(1, 1, 1, 1,
                                               1, 1), (totdim - 3)/6 - 1))/sqrt(totdim)
    ksi2 <- c(-1/2, -1/2, 1, -1/2, -1/2, 1, -1/2, -1/2, 1,
              rep(c(-1/2, -1/2, 1, -1/2, -1/2, 1), (totdim - 3)/6 -
                    1))/sqrt(3 * ((totdim - 3)/6 - 1) + 9/2)
    ksi3 <- c(-1, 1, 0, -1, 1, 0, -1, 1, 0, rep(c(-1, 1,
                                                  0, -1, 1, 0), (totdim - 3)/6 - 1))/sqrt(4 * ((totdim -
                                                                                                  3)/6 - 1) + 6)
    ksi4 <- c(-1/2, -1/2, -1/2, -1/2, -1/2, -1/2, 1, 1, 1,
              rep(c(1, 1, 1, -1, -1, -1), (totdim - 3)/6 - 1))/sqrt(6 *
                                                                      ((totdim - 3)/6 - 1) + 9/2)
    ksi5 <- c(1/4, 1/4, -1/2, 1/4, 1/4, -1/2, -1/2, -1/2,
              1, rep(c(-1/2, -1/2, 1, 1/2, 1/2, -1), (totdim -
                                                        3)/6 - 1))/sqrt(3 * ((totdim - 3)/6 - 1) + 9/4)
    ksi6 <- c(1/2, -1/2, 0, 1/2, -1/2, 0, -1, 1, 0, rep(c(-1,
                                                          1, 0, 1, -1, 0), (totdim - 3)/6 - 1))/sqrt(4 * ((totdim -
                                                                                                             3)/6 - 1) + 3)
    ksi <- cbind(ksi1, ksi2, ksi3, ksi4, ksi5, ksi6)[1:totdim,
                                                     1:(ncomp+1)]
  }
  if ((totdim%%6 == 4) & (totdim >= 6)) {
    dimok <- TRUE
    if (!(ncomp %in% 1:6)) {
      warning(paste("ncomp must be <= ", 6, "\n"))
      warning(paste("ncomp was set to ", 6, "\n"))
      ncomp <- 6
    }
    ksi1 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, rep(c(1, 1, 1,
                                                  1, 1, 1), (totdim - 4)/6 - 1))/sqrt(totdim)
    ksi2 <- c(1, -1, 1, -1, 0, 1, -1, 1, -1, 0, rep(c(-1/2,
                                                      -1/2, 1, -1/2, -1/2, 1), (totdim - 4)/6 - 1))/sqrt(3 *
                                                                                                           ((totdim - 4)/6 - 1) + 8)
    ksi3 <- c(1, 1, -1, -1, 0, 1, 1, -1, -1, 0, rep(c(-1,
                                                      1, 0, -1, 1, 0), (totdim - 4)/6 - 1))/sqrt(4 * ((totdim -
                                                                                                         4)/6 - 1) + 8)
    ksi4 <- c(1, -1, -1, 1, 0, 1, -1, -1, 1, 0, rep(c(1,
                                                      1, 1, -1, -1, -1), (totdim - 4)/6 - 1))/sqrt(6 *
                                                                                                     ((totdim - 4)/6 - 1) + 8)
    ksi5 <- c(1/4, 1/4, 1/4, 1/4, -1, 1/4, 1/4, 1/4, 1/4,
              -1, rep(c(-1/2, -1/2, 1, 1/2, 1/2, -1), (totdim -
                                                         4)/6 - 1))/sqrt(3 * ((totdim - 4)/6 - 1) + 5/2)
    ksi6 <- c(1, 1, 1, 1, 1, -1, -1, -1, -1, -1, rep(c(-1,
                                                       1, 0, 1, -1, 0), (totdim - 4)/6 - 1))/sqrt(4 * ((totdim -
                                                                                                          4)/6 - 1) + 10)
    ksi <- cbind(ksi1, ksi2, ksi3, ksi4, ksi5, ksi6)[1:totdim,
                                                     1:(ncomp+1)]
  }
  if ((totdim%%6 == 5) & (totdim >= 6)) {
    dimok <- TRUE
    if (!(ncomp %in% 1:6)) {
      warning(paste("ncomp must be <= ", 6, "\n"))
      warning(paste("ncomp was set to ", 6, "\n"))
      ncomp <- 6
    }
    ksi1 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, rep(c(1, 1,
                                                     1, 1, 1, 1), (totdim - 5)/6 - 1))/sqrt(totdim)
    ksi2 <- c(1, -1, 1, -1, 0, 1, -1, 1, -1, 0, 0, rep(c(-1/2,
                                                         -1/2, 1, -1/2, -1/2, 1), (totdim - 5)/6 - 1))/sqrt(3 *
                                                                                                              ((totdim - 5)/6 - 1) + 8)
    ksi3 <- c(1, 1, -1, -1, 0, 1, 1, -1, -1, 0, 0, rep(c(-1,
                                                         1, 0, -1, 1, 0), (totdim - 5)/6 - 1))/sqrt(4 * ((totdim -
                                                                                                            5)/6 - 1) + 8)
    ksi4 <- c(10/11, -12/11, -12/11, 10/11, -1/11, 10/11,
              -12/11, -12/11, 10/11, -1/11, 10/11, rep(c(1, 1,
                                                         1, -1, -1, -1), (totdim - 5)/6 - 1))/sqrt(6 *
                                                                                                     ((totdim - 5)/6 - 1) + 98/11)
    ksi5 <- c(13/196, 53/196, 53/196, 13/196, -53/49, 13/196,
              53/196, 53/196, 13/196, -53/49, 40/49, rep(c(-1/2,
                                                           -1/2, 1, 1/2, 1/2, -1), (totdim - 5)/6 - 1))/sqrt(3 *
                                                                                                               ((totdim - 5)/6 - 1) + 325/98)
    ksi6 <- c(4/5, 62/65, 62/65, 4/5, 77/65, -6/5, -68/65,
              -68/65, -6/5, -53/65, 8/13, rep(c(-1, 1, 0, 1, -1,
                                                0), (totdim - 5)/6 - 1))/sqrt(4 * ((totdim -
                                                                                      5)/6 - 1) + 138/13)
    ksi <- cbind(ksi1, ksi2, ksi3, ksi4, ksi5, ksi6)[1:totdim,
                                                     1:(ncomp+1)]
  }
  if (!dimok) {
    stop("Incorrect value for 'totdim'. 'totdim' must be > 1")
  }
  epsilon <- stats::rnorm(totdim, mean = rep(0, totdim), sd = varepsilon)
  #r <- stats::rnorm(6, mean = rep(0, 6), sd = varsR[1:6])
  r <- stats::rexp(6,varsR[1:6])
  simX <- r[1:(ncomp+1)] %*% t(ksi) + epsilon
  if ((ncomp+1) == 2) {
    HH <- 3
    eta21 <- c(1, 1, 1)/sqrt(3)
    eta22 <- c(1, 1, 1)/sqrt(3)
    eta <- cbind(eta21, eta22)
  }
  if ((ncomp+1) == 3) {
    HH <- 3
    eta31 <- c(1, 1, 1)/sqrt(3)
    eta32 <- c(1, 1, 1)/sqrt(3)
    eta33 <- c(1, 1, 1)/sqrt(3)
    eta <- cbind(eta31, eta32, eta33)
  }
  if ((ncomp+1) == 4) {
    HH <- 4
    eta41 <- c(1, 1, 1, 1)/2
    eta42 <- c(1, 1, 1, 1)/2
    eta43 <- c(1, 1, 1, 1)/2
    eta44 <- c(1, 1, 1, 1)/2
    eta <- cbind(eta41, eta42, eta43, eta44)
  }
  if ((ncomp+1) == 5) {
    HH <- 4
    eta51 <- c(1, 1, 1, 1)/2
    eta52 <- c(1, 1, 1, 1)/2
    eta53 <- c(1, 1, 1, 1)/2
    eta54 <- c(1, 1, 1, 1)/2
    eta55 <- c(1, 1, 1, 1)/2
    eta <- cbind(eta51, eta52, eta53, eta54, eta55)
  }
  if ((ncomp+1) == 6) {
    HH <- 4
    eta61 <- c(1, 1, 1, 1)/2
    eta62 <- c(1, 1, 1, 1)/2
    eta63 <- c(1, 1, 1, 1)/2
    eta64 <- c(1, 1, 1, 1)/2
    eta65 <- c(1, 1, 1, 1)/2
    eta66 <- c(1, 1, 1, 1)/2
    eta <- cbind(eta61, eta62, eta63, eta64, eta65, eta66)
  }
  #f <- stats::rnorm(6, mean = rep(0, 6), sd = varsF[1:6])
  f <- stats::rexp(6, varsF[1:6])
  z <- f[c(1:ncomp,(ncomp+2))] + r[c(1:ncomp,(ncomp+2))]
  sigmaScarre <- 0.001
  lambda <- 0.6
  sigmaPsi <- sigmaScarre * ((1 - lambda) * diag(rep(1, HH)) +
                               lambda * rep(1, HH) %*% t(rep(1, HH)))
  Psi <- mvtnorm::rmvnorm(1, mean = rep(0, HH), sigma = sigmaPsi)
  if ((z %*% t(eta) + Psi + offset)[1]<0){
    linearYgamma <- -((z %*% t(eta) + Psi + offset)[1])
    #simX=rep(-1,totdim)*simX
  }
  else {
    linearYgamma <- (z %*% t(eta) + Psi + offset)[1]
  }
  pYgamma <- linkobj$linkinv(linearYgamma)
  Ygamma <- rgamma(1, shape = pYgamma*10, scale=1/10)
  res <- c(Ygamma, simX)
  names(res) <- c("Ygamma", paste("X", 1:totdim, sep = ""))
  return(res)
}
