#' Draw samples from a zero-inflated poisson distribution
#'
#' @param n the number of samples to draw
#' @param lambda The poisson rate parameter
#' @param pstr0 probability of drawing a zero
#' @return Poisson counts of length \eqn{n}
#' @importFrom stats qpois dpois runif
#' @export
rzipois <- function(n, lambda, pstr0 = 0) {
    ans <- rpois(n, lambda)
    ans <- ifelse(runif(n) < pstr0, 0, ans)
    prob0 <- exp(-lambda)
    deflat.limit <- -1/expm1(lambda)
    ind0 <- (deflat.limit <= pstr0) & (pstr0 < 0)
    if (any(ind0)) {
        pobs0 <- pstr0[ind0] + (1 - pstr0[ind0]) * pstr0[ind0]
        ans[ind0] <- qpois(p = runif(sum(ind0),
                    min = dpois(0, lambda[ind0])), lambda[ind0])
        ans[ind0] <- ifelse(runif(sum(ind0)) < pobs0, 0, ans[ind0])
    }
    ans[pstr0 < deflat.limit] <- NaN
    ans[pstr0 > 1] <- NaN
    ans
}


#' @keywords internal
.zipois_getLam <- function(mu, S) {
    S <- max(sqrt(mu), S)
    (S^2/mu) + mu - 1
}

#' @keywords internal
.zipois_getP <- function(mu, S) {
    S <- max(sqrt(mu), S)
    (S^2 - mu) / (mu^2 - mu + S^2)
}

#' Generate multivariate, Zero-inflated poisson data,
#' with counts approximately correlated according to Sigma
#'
#' @param n number of samples to draw
#' @param mu mean vector for variables (of length \eqn{D})
#' @param Sigma \eqn{DxD} covariance or correlation matrix
#' @param lambdas supply rate parameter (instead of mu)
#' @param ps probability of zeros (instead of mu)
#' @param ... arguments passed to \code{VGAM::qzipois}
#' @return \eqn{Dxn} matrix with zi-poisson data
#' @importFrom VGAM qzipois
#' @export
rmvzipois <- function(n, mu, Sigma=diag(length(mu)), lambdas, ps, ...) {
    d   <- ncol(Sigma)
    Cor <- cov2cor(Sigma)
    SDs <- sqrt(diag(Sigma))

    if (missing(lambdas) || missing(ps)) {
        if (missing(mu)) stop("Need to supply mu")
        if (length(mu) != length(SDs)) stop("Sigma and mu dimensions don't match")
        lambdas <- unlist(lapply(1:length(SDs), function(i) .zipois_getLam(mu[i], SDs[i])))
        ps   <- unlist(lapply(1:length(SDs), function(i) .zipois_getP(mu[i], SDs[i])))
    }
    if (length(lambdas) != length(SDs)) stop("Sigma and mu/lambdas dimensions don't match")
    if (length(lambdas) == 1) stop("Need more than 1 variable")

    normd  <- rmvnorm(n, rep(0, d), Cor)
    unif   <- pnorm(normd)
    data <- matrix(VGAM::qzipois(unif, lambdas, pstr0=ps, ...), n, d)
    data <- .fixInf(data)
    return(data)
}


#' Generate multivariate poisson data,
#' with counts approximately correlated according to Sigma
#'
#' @param n number of samples to draw
#' @param mu mean vector for variables (of length \eqn{D})
#' @param Sigma \eqn{DxD} covariance or correlation matrix
#' @param ... Arguments passed to \code{qpois}
#' @return \eqn{Dxn} matrix with zi-poisson data
#' @importFrom stats qpois
#' @export
rmvpois <- function(n, mu, Sigma, ...) {
    Cor <- cov2cor(Sigma)
    SDs <- sqrt(diag(Sigma))
    d   <- length(SDs)
    if (length(mu) != length(SDs)) stop("Sigma and mu/lambdas dimensions don't match")
    if (length(mu) == 1) stop("Need more than 1 variable")
    normd  <- rmvnorm(n, rep(0, d), Cor)
    unif   <- pnorm(normd)
    data <- matrix(qpois(unif, mu, ...), n, d)
    data <- .fixInf(data)
    return(data)
}

#' @keywords internal
.negbin_getK <- function(mu, S) {
    S <- max(sqrt(mu), S)
    mu^2/((S^2+1e-3)-mu)
}

#' Generate multivariate, Zero-inflated negative binomial data,
#' with counts approximately correlated according to Sigma
#'
#' @param n number of samples to draw
#' @param mu mean vector for variables (of length \eqn{D})
#' @param Sigma \eqn{DxD} covariance or correlation matrix
#' @param ks shape parameter
#' @param ... other arguments to the negative binomial distribution
#' @return \eqn{Dxn} matrix with zi-poisson data
#' @importFrom stats qnbinom
#' @export
rmvnegbin <- function(n, mu, Sigma, ks, ...) {
# Generate an NxD matrix of Zero-inflated poisson data,
# with counts approximately correlated according to Sigma
    Cor <- cov2cor(Sigma)
    SDs <- sqrt(diag(Sigma))
    if (missing(mu)) stop('mu is required')
    if (length(mu) != length(SDs)) stop("Sigma and mu dimensions don't match")
    if (missing(ks)) {
        ks   <- unlist(lapply(1:length(SDs), function(i) .negbin_getK(mu[i], SDs[i])))
    }
    d   <- length(mu)
    normd  <- rmvnorm(n, rep(0, d), Sigma=Cor)
    unif   <- pnorm(normd)
    data <- t(qnbinom(t(unif), mu=mu, size=ks, ...))
    data <- .fixInf(data)
    return(data)
}


#' @keywords internal
.zinegbin_getLam <- function(mu, S) {
    S   <- max(sqrt(mu)+1e-3, S)
    (mu + (mu^2 - mu + S^2) / mu) / 2
}

#' @keywords internal
.zinegbin_getP <- function(mu, lam) {
    1 - (mu / lam)
}

#' @keywords internal
.zinegbin_getK <- function(mu, S, lam) {
    S   <- max(sqrt(mu)+1e-3, S)
    (mu * lam) / (mu^2 - (mu * (lam + 1)) + S^2)
}


#' Generate multivariate, negative binomial data,
#' with counts approximately correlated according to Sigma
#'
#' @param n number of samples to draw
#' @param mu mean vector for variables (of length \eqn{D})
#' @param Sigma \eqn{DxD} covariance or correlation matrix
#' @param ps probability of zero inflation
#' @param munbs Rate/mean parameter (instead of mu)
#' @param ks shape parameter
#' @param ... other arguments to the negative binomial distribution
#' @return \eqn{Dxn} matrix with zi-poisson data
#' @importFrom VGAM qzinegbin
#' @export
rmvzinegbin <- function(n, mu, Sigma, munbs, ks, ps, ...) {
# Generate an NxD matrix of Zero-inflated poisson data,
# with counts approximately correlated according to Sigma
    Cor <- cov2cor(Sigma)
    SDs <- sqrt(diag(Sigma))
    if (missing(munbs) || missing(ps) || missing(ks)) {
        if (length(mu) != length(SDs)) stop("Sigma and mu dimensions don't match")
        munbs <- unlist(lapply(1:length(SDs), function(i) .zinegbin_getLam(mu[i], SDs[i])))
        ps   <- unlist(lapply(1:length(SDs), function(i) .zinegbin_getP(mu[i], munbs[i])))
        ks   <- unlist(lapply(1:length(SDs), function(i) .zinegbin_getK(mu[i], SDs[i], munbs[i])))
    }
    if (length(munbs) != length(SDs)) stop("Sigma and mu dimensions don't match")
    d   <- length(munbs)
    normd  <- rmvnorm(n, rep(0, d), Sigma=Cor)
    unif   <- pnorm(normd)
    data <- matrix(VGAM::qzinegbin(unif, munb=munbs, size=ks, pstr0=ps, ...), n, d)
    data <- .fixInf(data)
    return(data)
}



#' @keywords internal
.fixInf <- function(data) {
    # hacky way of replacing infinite values with the col max + 1
    if (any(is.infinite(data))) {
       data <-  apply(data, 2, function(x) {
              if (any(is.infinite(x))) {
                   x[ind<-which(is.infinite(x))] <- NA
                   x[ind] <- max(x, na.rm=TRUE)+1
                 }
                x
                })
    }
    data
}


#' Draw samples from multivariate, correlated normal distribution
#' with counts correlated according to Sigma
#'
#' @param n number of samples to draw
#' @param mu mean vector for variables (of length \eqn{D})
#' @param Sigma \eqn{DxD} covariance or correlation matrix
#' @param tol numerical tolerance for a zero eigenvalue (check for PD of Sigma)
#' @param empirical is Sigma the empirical correlation?
#' @return \eqn{Dxn} matrix with Gaussian data
#' @export
rmvnorm <- function(n=100, mu=rep(0,10), Sigma=diag(10), tol=1e-6, empirical=TRUE) {
    p <- length(mu)
    if (!all(dim(Sigma) == c(p, p)))
        stop("incompatible arguments")
    eS <- eigen(Sigma, symmetric = TRUE)
    ev <- eS$values
    if (!all(ev >= -tol * abs(ev[1L])))
        stop("'Sigma' is not positive definite")
    X <- matrix(rnorm(p * n), n)
    if (empirical) {
        X <- scale(X, TRUE, FALSE)
        X <- X %*% svd(X, nu = 0, nv = length(mu))$v
        X <- scale(X, FALSE, TRUE)
    }
    X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
    return(t(X))
}

#' Convert a symmetric correlation matrix to a covariance matrix
#' given the standard deviation
#'
#' @param cor a symmetric correlation matrix
#' @param sds standard deviations of the resulting covariance.
#' @return Covariance matrix of sample dimension as cor
#' @export
cor2cov <- function(cor, sds) {
    if (length(sds) != length(diag(cor))) stop("inputs are of mismatched dimension")
    cor * sds * rep(sds, each=nrow(cor))
}
