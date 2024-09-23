#' Non-parametric (Y,T) Bootstrap for selecting the number of components in PLS 
#' GLR models
#'
#' Provides a wrapper for the bootstrap function \code{boot} from the
#' \code{boot} R package.\cr Implements non-parametric bootstraps for PLS
#' Generalized Linear Regression models by (Y,T) resampling to select the 
#' number of components.
#' 
#' More details on bootstrap techniques are available in the help of the
#' \code{\link[boot:boot]{boot}} function.
#' 
#' @param object An object of class \code{plsRmodel} to bootstrap
#' @param typeboot The type of bootstrap. (\code{typeboot="boot_comp"}) for 
#' (Y,T) bootstrap to select components. Defaults to 
#' (\code{typeboot="boot_comp"}).
#' @param R The number of bootstrap replicates. Usually this will be a single
#' positive integer. For importance resampling, some resamples may use one set
#' of weights and others use a different set of weights. In this case \code{R}
#' would be a vector of integers where each component gives the number of
#' resamples from each of the rows of weights.
#' @param statistic A function which when applied to data returns a vector
#' containing the statistic(s) of interest. \code{statistic} must take at least
#' two arguments. The first argument passed will always be the original data.
#' The second will be a vector of indices, frequencies or weights which define
#' the bootstrap sample. Further, if predictions are required, then a third
#' argument is required which would be a vector of the random indices used to
#' generate the bootstrap predictions. Any further arguments can be passed to
#' statistic through the \code{...} argument.
#' @param sim A character string indicating the type of simulation required.
#' Possible values are \code{"ordinary"} (the default), \code{"balanced"},
#' \code{"permutation"}, or \code{"antithetic"}.
#' @param stype A character string indicating what the second argument of
#' \code{statistic} represents. Possible values of stype are \code{"i"}
#' (indices - the default), \code{"f"} (frequencies), or \code{"w"} (weights).
#' @param stabvalue A value to hard threshold bootstrap estimates computed from
#' atypical resamplings. Especially useful for Generalized Linear Models.
#' @param \dots Other named arguments for \code{statistic} which are passed
#' unchanged each time it is called. Any such arguments to \code{statistic}
#' should follow the arguments which \code{statistic} is required to have for
#' the simulation. Beware of partial matching to arguments of \code{boot}
#' listed above.
#' @return An object of class \code{"boot"}. See the Value part of the help of
#' the function \code{\link[boot:boot]{boot}}.
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
#' \doi{10.3389/fams.2021.693126}\cr.
#' 
#' @examples
#' set.seed(314)
#' library(plsRglm)
#' data(aze_compl, package="plsRglm")
#' Xaze_compl<-aze_compl[,2:34]
#' yaze_compl<-aze_compl$y
#' dataset <- cbind(y=yaze_compl,Xaze_compl)
#' modplsglm <- plsRglm::plsRglm(y~.,data=dataset,10,modele="pls-glm-family", family = binomial)
#' 
#' comp_aze_compl.bootYT <- nbcomp.bootplsRglm(modplsglm, R=250)
#' boxplots.bootpls(comp_aze_compl.bootYT)
#' confints.bootpls(comp_aze_compl.bootYT)
#' plots.confints.bootpls(confints.bootpls(comp_aze_compl.bootYT),typeIC = "BCa")
#' 
#' comp_aze_compl.permYT <- nbcomp.bootplsRglm(modplsglm, R=250, sim="permutation")
#' boxplots.bootpls(comp_aze_compl.permYT)
#' confints.bootpls(comp_aze_compl.permYT, typeBCa=FALSE)
#' plots.confints.bootpls(confints.bootpls(comp_aze_compl.permYT, typeBCa=FALSE))
nbcomp.bootplsRglm <- function (object, typeboot = "boot_comp", R = 250, statistic = coefs.plsRglm.CSim,
                                sim = "ordinary", stype = "i", stabvalue = 1e+06, ...)
{
  callplsRglm <- object$call
  maxcoefvalues <- stabvalue * abs(object$Coeffs)
  dataset <- cbind(y = object$dataY, object$dataX)
  nt <- eval(callplsRglm$nt)
  ifbootfail <- as.matrix(as.numeric(ifelse(any(class(dataset[,
                                                              1]) == "factor"), rep(NA, ncol(dataset) + nlevels(dataset[,
                                                                                                                        1]) - 1), rep(NA, ncol(dataset)))))
  if (!is.null(callplsRglm$modele)) {
    modele <- eval(callplsRglm$modele)
  }
  else {
    modele <- "pls"
  }
  if (!is.null(callplsRglm$family)) {
    family <- eval(callplsRglm$family)
  }
  else {
    family <- NULL
  }
  # if (typeboot == "plsmodel") {
  #   temp.bootplsRglm <- if (!(sim == "permutation")) {
  #     boot::boot(data = dataset, statistic = coefs.plsRglm, sim = sim,
  #          stype = stype, R = R, nt = nt, modele = modele,
  #          family = family, maxcoefvalues = maxcoefvalues,
  #          ifbootfail = ifbootfail, ...)
  #   }
  #   else {
  #     boot::boot(data = dataset, statistic = permcoefs.plsRglm,
  #          sim = sim, stype = stype, R = R, nt = nt, modele = modele,
  #          family = family, maxcoefvalues = maxcoefvalues,
  #          ifbootfail = ifbootfail)
  #   }
  #   indices.temp.bootplsRglm <- !is.na(temp.bootplsRglm$t[,
  #                                                         1])
  #   temp.bootplsRglm$t = temp.bootplsRglm$t[indices.temp.bootplsRglm,
  #   ]
  #   temp.bootplsRglm$R = sum(indices.temp.bootplsRglm)
  #   temp.bootplsRglm$call$R <- sum(indices.temp.bootplsRglm)
  #   return(temp.bootplsRglm)
  # }
  # if (typeboot == "fmodel_np") {
  #   dataRepYtt <- cbind(y = object$RepY, object$tt)
  #   wwetoile <- object$wwetoile
  #   temp.bootplsRglm <- if (!(sim == "permutation")) {
  #     boot::boot(data = dataRepYtt, statistic = coefs.plsRglmnp,
  #          sim = sim, stype = stype, R = R, nt = nt, modele = modele,
  #          family = family, maxcoefvalues = maxcoefvalues[-(1:(length(object$Coeffs) -
  #                                                                ncol(object$dataX)))], wwetoile = wwetoile,
  #          ifbootfail = ifbootfail, ...)
  #   }
  #   else {
  #     boot::boot(data = dataRepYtt, statistic = permcoefs.plsRglmnp,
  #          sim = sim, stype = stype, R = R, nt = nt, modele = modele,
  #          family = family, maxcoefvalues = maxcoefvalues[-(1:(length(object$Coeffs) -
  #                                                                ncol(object$dataX)))], wwetoile = wwetoile,
  #          ifbootfail = ifbootfail)
  #   }
  #   indices.temp.bootplsRglm <- !is.na(temp.bootplsRglm$t[,
  #                                                         1])
  #   temp.bootplsRglm$t = temp.bootplsRglm$t[indices.temp.bootplsRglm,
  #   ]
  #   temp.bootplsRglm$R = sum(indices.temp.bootplsRglm)
  #   temp.bootplsRglm$call$R <- sum(indices.temp.bootplsRglm)
  #   return(temp.bootplsRglm)
  # }
  if (typeboot == "boot_comp") {
    maxcoefvalues <- stabvalue * abs(object$CoeffC)
    dataRepYtt <- cbind(y = object$RepY, object$tt)
    temp.bootplsRglm <- if (!(sim == "permutation")) {
      boot::boot(data = dataRepYtt, statistic = coefs.plsRglm.CSim,
                 sim = sim, stype = stype, R = R, nt = nt, modele = modele,
                 family = family, maxcoefvalues = maxcoefvalues,
                 ifbootfail = ifbootfail, ...)
    }
    else {
      boot::boot(data = dataRepYtt, statistic = permcoefs.plsRglm.CSim,
                 sim = sim, stype = stype, R = R, nt = nt, modele = modele,
                 family = family, maxcoefvalues = maxcoefvalues[-(1:(length(object$Coeffs) -
                                                                       ncol(object$dataX)))],
                 ifbootfail = ifbootfail)
    }
    indices.temp.bootplsRglm <- !is.na(temp.bootplsRglm$t[,
                                                          1])
    temp.bootplsRglm$t = temp.bootplsRglm$t[indices.temp.bootplsRglm,
    ]
    temp.bootplsRglm$R = sum(indices.temp.bootplsRglm)
    temp.bootplsRglm$call$R <- sum(indices.temp.bootplsRglm)
    return(temp.bootplsRglm)
  }
  # if (typeboot == "fmodel_par") {
  #   temp.bootplsRglm <- if (!(sim == "permutation")) {
  #     boot::boot(data = dataset, statistic = coefs.plsRglm, sim = sim,
  #          stype = stype, R = R, nt = nt, modele = modele,
  #          family = family, maxcoefvalues = maxcoefvalues,
  #          ifbootfail = ifbootfail, ...)
  #   }
  #   else {
  #     boot::boot(data = dataset, statistic = permcoefs.plsRglm,
  #          sim = sim, stype = stype, R = R, nt = nt, modele = modele,
  #          family = family, maxcoefvalues = maxcoefvalues,
  #          ifbootfail = ifbootfail)
  #   }
  #   indices.temp.bootplsRglm <- !is.na(temp.bootplsRglm$t[,
  #                                                         1])
  #   temp.bootplsRglm$t = temp.bootplsRglm$t[indices.temp.bootplsRglm,
  #   ]
  #   temp.bootplsRglm$R = sum(indices.temp.bootplsRglm)
  #   temp.bootplsRglm$call$R <- sum(indices.temp.bootplsRglm)
  #   return(temp.bootplsRglm)
  # }
}
