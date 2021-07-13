#' @title bootPLS-package
#'
#' @description: Several implementations of non-parametric stable bootstrap-based techniques to determine the numbers of components for Partial Least Squares linear or generalized linear regression models as well as and sparse Partial Least Squares linear or generalized linear regression models. The package collect techniques that were published in a book chapter (Magnanensi et al. 2016, 'The Multiple Facets of Partial Least Squares and Related Methods', <\doi{10.1007/978-3-319-40643-5_18}) and two articles (Magnanensi et al. 2017, 'Statistics and Computing', <\doi{10.1007/s11222-016-9651-4}) and (Magnanensi et al. 2021, 'Frontiers in Applied Mathematics and Statistics', \doi{10.3389/fams.2021.693126}).
#'
#' @docType package
#' @name bootPLS-package
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
#' @importFrom graphics par strwidth
#' @importFrom stats binomial coef glm hat lm make.link median predict rgamma sd uniroot weighted.mean
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom foreach foreach %dopar%
#'            
#' @examples
#' set.seed(314)
#' library(bootPLS)
#' data(datasim)
#' head(datasim)
#' 
NULL
