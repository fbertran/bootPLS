#' @title Internal bigPLS functions
#' 
#' @name internal-bootPLS
#' 
#' @description These are not to be called by the user.
#' 
#' @aliases ust spls.dv correctp correctp.withoutK spls.Cboot cv.split
#' @author Jérémy Magnanensi, Frédéric Bertrand\cr
#' \email{frederic.bertrand@@utt.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
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
#' accepted.
#' 
#' @keywords internal
NULL


### For spls
ust<-function (b, eta)
{
  b.ust <- matrix(0, length(b), 1)
  if (eta < 1) {
    valb <- abs(b) - eta * max(abs(b))
    b.ust[valb >= 0] <- valb[valb >= 0] * (sign(b))[valb >=
                                                      0]
  }
  return(b.ust)
}

spls.dv<-function (Z, eta, kappa, eps, maxstep)
{
  p <- nrow(Z)
  q <- ncol(Z)
  Znorm1 <- median(abs(Z))
  Z <- Z/Znorm1
  if (q == 1) {
    c <- ust(Z, eta)
  }
  if (q > 1) {
    M <- Z %*% t(Z)
    dis <- 10
    i <- 1
    if (kappa == 0.5) {
      c <- matrix(10, p, 1)
      c.old <- c
      while (dis > eps & i <= maxstep) {
        mcsvd <- svd(M %*% c)
        a <- mcsvd$u %*% t(mcsvd$v)
        c <- ust(M %*% a, eta)
        dis <- max(abs(c - c.old))
        c.old <- c
        i <- i + 1
      }
    }
    if (kappa > 0 & kappa < 0.5) {
      kappa2 <- (1 - kappa)/(1 - 2 * kappa)
      c <- matrix(10, p, 1)
      c.old <- c
      h <- function(lambda) {
        alpha <- solve(M + lambda * diag(p)) %*% M %*%
          c
        obj <- t(alpha) %*% alpha - 1/kappa2^2
        return(obj)
      }
      if (h(eps) * h(1e+30) > 0) {
        while (h(eps) <= 1e+05) {
          M <- 2 * M
          c <- 2 * c
        }
      }
      while (dis > eps & i <= maxstep) {
        if (h(eps) * h(1e+30) > 0) {
          while (h(eps) <= 1e+05) {
            M <- 2 * M
            c <- 2 * c
          }
        }
        lambdas <- uniroot(h, c(eps, 1e+30))$root
        a <- kappa2 * solve(M + lambdas * diag(p)) %*%
          M %*% c
        c <- ust(M %*% a, eta)
        dis <- max(abs(c - c.old))
        c.old <- c
        i <- i + 1
      }
    }
  }
  return(c)
}

correctp=function (x, y, eta, K, kappa, select, fit) 
{
  if (min(eta) < 0 | max(eta) >= 1) {
    if (max(eta) == 1) {
      stop("eta should be strictly less than 1!")
    }
    if (length(eta) == 1) {
      stop("eta should be between 0 and 1!")
    }
    else {
      stop("eta should be between 0 and 1! \n  Choose appropriate range of eta!")
    }
  }
  if (max(K) > ncol(x)) {
    stop("K cannot exceed the number of predictors! Pick up smaller K!")
  }
  if (max(K) >= nrow(x)) {
    stop("K cannot exceed the sample size! Pick up smaller K!")
  }
  if (min(K) <= 0 | !all(K%%1 == 0)) {
    if (length(K) == 1) {
      stop("K should be a positive integer!")
    }
    else {
      stop("K should be a positive integer! \n  Choose appropriate range of K!")
    }
  }
  if (kappa > 0.5 | kappa < 0) {
    cat("kappa should be between 0 and 0.5! kappa=0.5 is used. \n\n")
    kappa <- 0.5
  }
  if (select != "pls2" & select != "simpls") {
    cat("Invalid PLS algorithm for variable selection.\n")
    cat("pls2 algorithm is used. \n\n")
    select <- "pls2"
  }
  fits <- c("simpls", "kernelpls", "widekernelpls", "oscorespls")
  if (!any(fit == fits)) {
    cat("Invalid PLS algorithm for model fitting\n")
    cat("simpls algorithm is used. \n\n")
    fit <- "simpls"
  }
  list(K = K, eta = eta, kappa = kappa, select = select, fit = fit)
}

correctp.withoutK=function (x, y, eta, kappa, select, fit) 
{
  if (min(eta) < 0 | max(eta) >= 1) {
    if (max(eta) == 1) {
      stop("eta should be strictly less than 1!")
    }
    if (length(eta) == 1) {
      stop("eta should be between 0 and 1!")
    }
    else {
      stop("eta should be between 0 and 1! \n  Choose appropriate range of eta!")
    }
  }
  if (kappa > 0.5 | kappa < 0) {
    cat("kappa should be between 0 and 0.5! kappa=0.5 is used. \n\n")
    kappa <- 0.5
  }
  if (select != "pls2" & select != "simpls") {
    cat("Invalid PLS algorithm for variable selection.\n")
    cat("pls2 algorithm is used. \n\n")
    select <- "pls2"
  }
  fits <- c("simpls", "kernelpls", "widekernelpls", "oscorespls")
  if (!any(fit == fits)) {
    cat("Invalid PLS algorithm for model fitting\n")
    cat("simpls algorithm is used. \n\n")
    fit <- "simpls"
  }
  list(eta = eta, kappa = kappa, select = select, fit = fit)
}

spls.Cboot=function (x, y, K, eta, kappa = 0.5, select = "pls2", fit = "simpls", 
                     scale.x = TRUE, scale.y = FALSE, eps = 1e-04, maxstep = 100, 
                     trace = FALSE) 
{
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ip <- c(1:p)
  y <- as.matrix(y)
  q <- ncol(y)
  one <- matrix(1, 1, n)
  mu <- one %*% y/n
  y <- scale(y, drop(mu), FALSE)
  meanx <- drop(one %*% x)/n
  x <- scale(x, meanx, FALSE)
  if (scale.x) {
    normx <- sqrt(drop(one %*% (x^2))/(n - 1))
    if (any(normx < .Machine$double.eps)) {
      stop("Some of the columns of the predictor matrix have zero variance.")
    }
    x <- scale(x, FALSE, normx)
  }
  else {
    normx <- rep(1, p)
  }
  if (scale.y) {
    normy <- sqrt(drop(one %*% (y^2))/(n - 1))
    if (any(normy < .Machine$double.eps)) {
      stop("Some of the columns of the response matrix have zero variance.")
    }
    y <- scale(y, FALSE, normy)
  }
  else {
    normy <- rep(1, q)
  }
  betahat <- matrix(0, p, q)
  betamat <- list()
  x1 <- x
  y1 <- y
  type <- correctp(x, y, eta, K, kappa, select, fit)
  eta <- type$eta
  K <- type$K
  kappa <- type$kappa
  select <- type$select
  fit <- type$fit
  if (is.null(colnames(x))) {
    xnames <- c(1:p)
  }
  else {
    xnames <- colnames(x)
  }
  new2As <- list()
  if (trace) {
    cat("The variables that join the set of selected variables at each step:\n")
  }
  for (k in 1:K) {
    Z <- t(x1) %*% y1
    what <- spls.dv(Z, eta, kappa, eps, maxstep)
    A <- unique(ip[what != 0 | betahat[, 1] != 0])
    new2A <- ip[what != 0 & betahat[, 1] == 0]
    xA <- x[, A, drop = FALSE]
    plsfit <- pls::plsr(y ~ xA, ncomp = min(k, length(A)), 
                        method = fit, scale = FALSE)
    betahat <- matrix(0, p, q)
    betahat[A, ] <- matrix(coef(plsfit), length(A), q)
    betamat[[k]] <- betahat
    pj <- plsfit$projection
    if (select == "pls2") {
      y1 <- y - x %*% betahat
    }
    if (select == "simpls") {
      pw <- pj %*% solve(t(pj) %*% pj) %*% t(pj)
      x1 <- x
      x1[, A] <- x[, A, drop = FALSE] - x[, A, drop = FALSE] %*% 
        pw
    }
    new2As[[k]] <- new2A
    if (trace) {
      if (length(new2A) <= 10) {
        cat(paste("- ", k, "th step (K=", k, "):\n", 
                  sep = ""))
        cat(xnames[new2A])
        cat("\n")
      }
      else {
        cat(paste("- ", k, "th step (K=", k, "):\n", 
                  sep = ""))
        nlines <- ceiling(length(new2A)/10)
        for (i in 0:(nlines - 2)) {
          cat(xnames[new2A[(10 * i + 1):(10 * (i + 1))]])
          cat("\n")
        }
        cat(xnames[new2A[(10 * (nlines - 1) + 1):length(new2A)]])
        cat("\n")
      }
    }
  }
  coeffC <- pls::Yloadings(plsfit)[,1:min(K, length(A))]
  tt <- pls::scores(plsfit)[,1:min(K, length(A))]
  if (!is.null(colnames(x))) {
    rownames(betahat) <- colnames(x)
  }
  if (q > 1 & !is.null(colnames(y))) {
    colnames(betahat) <- colnames(y)
  }
  object <- list(x = x, y = y, coeffC = coeffC, tt = tt, betahat = betahat, A = A, betamat = betamat, 
                 new2As = new2As, mu = mu, meanx = meanx, normx = normx, 
                 normy = normy, eta = eta, K = K, kappa = kappa, select = select, 
                 fit = fit, projection = pj)
  class(object) <- "spls"
  object
}

cv.split=function (y, fold) 
{
  n <- length(y)
  group <- table(y)
  x <- c()
  for (i in 1:length(group)) {
    x.group <- c(1:n)[y == names(group)[i]]
    x <- c(x, sample(x.group))
  }
  foldi <- split(x, rep(1:fold, length = n))
  return(foldi)
}

wpls=function (x, y, V, K = ncol(x), type = "pls1", center.x = TRUE, 
               scale.x = FALSE) 
{
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  x1 <- x
  y1 <- y
  W <- matrix(0, p, K)
  T <- matrix(0, n, K)
  Q <- matrix(0, q, K)
  P <- matrix(0, p, K)
  for (k in 1:K) {
    w <- t(x1) %*% as.matrix(V * y1)
    w <- w/sqrt(sum(w^2))
    W[, k] <- w
    t <- x1 %*% w
    T[, k] <- t
    coef.q <- sum(t * V * y1)/sum(t * V * t)
    Q[, k] <- coef.q
    coef.p <- t(as.matrix(t * V)) %*% x1/sum(t * V * t)
    P[, k] <- coef.p
    if (type == "pls1") {
      y1 <- y1 - t %*% coef.q
      x1 <- x1 - t %*% coef.p
    }
    if (type == "simpls") {
      pj <- w
      pw <- pj %*% solve(t(pj) %*% pj) %*% t(pj)
      x1 <- x1 - x1 %*% pw
    }
  }
  list(W = W, T = T, Q = Q, P = P)
}

### Updating SGPLS function to get T

sgpls.T=function (x, y, K, eta, scale.x = TRUE, eps = 1e-05, denom.eps = 1e-20, 
                  zero.eps = 1e-05, maxstep = 100, br = TRUE, ftype = "iden") 
{
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  ip <- c(1:p)
  y <- as.matrix(y)
  q <- ncol(y)
  one <- matrix(1, 1, n)
  mu <- apply(x, 2, mean)
  x0 <- scale(x, mu, FALSE)
  if (scale.x) {
    sigma <- apply(x, 2, sd)
    x0 <- scale(x0, FALSE, sigma)
  }
  else {
    sigma <- rep(1, ncol(x))
    x0 <- x0
  }
  beta1hat <- matrix(0, p, q)
  beta1hat.old <- beta1hat + 1000
  beta0hat <- 0
  re <- 100
  min.re <- 1000
  nstep <- 0
  nstep.min <- 0
  while (re > eps & nstep < maxstep) {
    if (nstep == 0) {
      p0 <- (y + 0.5)/2
      V <- as.vector(p0 * (1 - p0))
      A <- c(1:p)
    }
    else {
      exp.xb <- exp(beta0hat + x0 %*% beta1hat)
      p0 <- exp.xb/(1 + exp.xb)
      p0[exp.xb == Inf] <- 1 - zero.eps
      p0[p0 < zero.eps] <- zero.eps
      p0[p0 > (1 - zero.eps)] <- 1 - zero.eps
      V <- as.vector(p0 * (1 - p0))
    }
    switch(ftype, hat = {
      H <- hat(sweep(cbind(rep(1, n), x0), 1, sqrt(V), 
                     "*"), intercept = FALSE)
    }, iden = {
      H <- rep(1, n)
    })
    if (nstep == 0) {
      y0 <- beta0hat + x0 %*% beta1hat + (y - p0)/V
    }
    else {
      V <- V * (H * br + 1)
      y0 <- beta0hat + x0 %*% beta1hat + (y + H * br/2 - 
                                            (H * br + 1) * p0)/V
    }
    y1 <- y0
    y1 <- y1 - mean(y1)
    x1 <- x0
    A.old <- c()
    for (k in 1:K) {
      Z <- t(x1) %*% as.matrix(V * y1)
      Znorm1 <- median(abs(Z))
      Z <- Z/Znorm1
      what <- ust(Z, eta)
      A <- sort(unique(c(A.old, ip[what != 0])))
      x0A <- x0[, A, drop = FALSE]
      plsfit <- wpls(x0A, y0, V, K = min(k, length(A)), 
                     type = "pls1", center.x = FALSE, scale.x = FALSE)
      A.old <- A
      y1 <- y0 - plsfit$T %*% t(plsfit$Q)
      x1 <- x0
      x1[, A] <- x0[, A] - plsfit$T %*% t(plsfit$P)
    }
    x0A <- x0[, A, drop = FALSE]
    plsfit <- wpls(x0A, y0, V, K = min(K, length(A)), type = "pls1", 
                   center.x = FALSE, scale.x = FALSE)
    W <- plsfit$W
    T <- plsfit$T
    P <- plsfit$P
    Q <- plsfit$Q
    beta1hat.old <- beta1hat
    beta1hat <- matrix(0, p, q)
    beta1hat[A, ] <- W %*% solve(t(P) %*% W) %*% t(Q)
    beta0hat <- weighted.mean((y0 - T %*% t(Q)), sqrt(V))
    re <- mean(abs(beta1hat - beta1hat.old))/mean(abs(beta1hat.old) + 
                                                    denom.eps)
    nstep <- nstep + 1
    if (re < min.re & nstep > 1) {
      min.re <- re
      nstep.min <- nstep
      beta1hat.min <- beta1hat
      beta0hat.min <- beta0hat
      A.min <- A
      W.min <- W
    }
  }
  if (re > eps) {
    if (nstep.min > 0) {
      converged <- FALSE
      beta1hat <- beta1hat.min
      beta0hat <- beta0hat.min
      A <- A.min
      W <- W.min
    }
  }
  betahat <- matrix(c(beta0hat, beta1hat))
  if (!is.null(colnames(x))) {
    rownames(betahat) <- 1:nrow(betahat)
    rownames(betahat)[1] <- "intercept"
    rownames(betahat)[2:nrow(betahat)] <- colnames(x)
  }
  else {
    rownames(betahat) <- c(0, paste("x", 1:p, sep = ""))
    rownames(betahat)[1] <- "intercept"
  }
  object <- list(x = x, y = y, x0 = x0, eta = eta, K = K, CoeffC=Q, tt=T, betahat = betahat, 
                 A = A, W = W, mu = mu, sigma = sigma)
  class(object) <- "sgpls"
  object
}

