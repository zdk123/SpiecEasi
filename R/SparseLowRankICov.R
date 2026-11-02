#' Sparse plus Low Rank inverse covariance
#'
#' Select an inverse covariance matrix that is a sparse plus low rank decomposition.
#'
#' @param data the n x p data matrix
#' @param npn flag to first fit nonparametric normal transform to the data
#' @param verbose flag to turn on verbose output
#' @param cor flag to use correlation matrix as the input (default: false - uses covariance)
#' @param ... arguments to override default algorithm settings (see details)
#' @export
#' @details
#' This is a wrapper function for sparse plus low rank iCov estimations performed by a custom ADMM algorithm.
#'
#' Therefore, arguments \code{...} should be named. Typically, these are for specifying a penalty parameter, \code{lambda}, or the number of penalties to use.
#' By default 10 pentalties are used, ranging logarithmically between \code{lambda.min.ratio}*MAX and MAX.
#' Max is the theoretical upper bound on \code{lambda} and us \code{max|S|}, the maximum absolute value in the data correlation matrix.
#' \code{lambda.min.ratio} is 1e-3 by default. Lower values of \code{lambda} require more memory/cpu time to compute, and sometimes huge will throw an error.
#'
#' The argument \code{nlambda} determines the number of penalties - somewhere between 10-100 is usually good, depending on how the values of empirical correlation are distributed.#' @export
#'
#' One of \code{beta} (penalty for the nuclear norm) or \code{r} (number of ranks) should be supplied or \code{r=2} is chosen by default.
#' @examples
#' # simulate data with 1 negative correlation
#'  set.seed(10010)
#'  Sigma <- diag(10)*2
#'  Sigma[1,2] <- Sigma[2,1] <- -.9
#'  data  <- exp(rmvnorm(50, runif(10, 0, 2), Sigma))
#'
#' # normalize
#'  data.f   <- t(apply(data, 1, norm_to_total))
#'  data.clr <- t(clr(data.f, 1))
#'
#' # estimate
#'  est.clr  <- sparseLowRankiCov(data.clr, cor=TRUE, r=2)
#'  est.f    <- sparseLowRankiCov(data.f, cor=TRUE, r=2)
#'  est.log  <- sparseLowRankiCov(log(data), cor=TRUE, r=2)
#'
#' # visualize results
#'  par(mfrow=c(1,3))
#'  image(as.matrix(est.log$path[[6]][1:5,1:5]))
#'  image(as.matrix(est.clr$path[[6]][1:5,1:5]))
#'  image(as.matrix(est.f$path[[6]][1:5,1:5]))
sparseLowRankiCov <- function(data, npn=FALSE, verbose=FALSE, cor=FALSE, ...) {
## TODO: make args to admm2 explicit
  args <- list(...)
  if (length(args$r) > 1 || length(args$beta) > 1)
    stop("Only single value allowed for \'r\' or \'beta\'")

  if (npn) data <- huge::huge.npn(data, verbose=verbose)
  if (isSymmetric(data)) SigmaO <- data
  else SigmaO <- cov(data)
  if (cor) SigmaO <- cov2cor(SigmaO)

  if (!is.null(args[[ "lambda.max" ]])) maxlam <- args$lambda.max
  else maxlam <- 1
  if (is.null(args[[ "lambda" ]])) {
    if (is.null(args[[ "lambda.min.ratio" ]])) args$lambda.min.ratio <- 1e-3
    if (is.null(args[[ "nlambda" ]])) args$nlambda <- 10
    lambda <- getLamPath(maxlam, maxlam*args$lambda.min.ratio, args$nlambda, log=TRUE)
    args$lambda.min.ratio <- NULL ; args$nlambda <- NULL
  } else lambda <- args$lambda
  args$lambda.min.ratio <- args$nlambda <- args$lambda <- args$lambda.max <- NULL

  if (is.null(args[[ "beta" ]])) {
    if (is.null(args[[ "r" ]]))
      args$r <- 2
  }

  n <- length(lambda)
  args$SigmaO <- SigmaO
###  args$r <- Lrank

#  lest <- parallel::mclapply(lambda, mc.cores=ncores, FUN=function(lam)
#                             do.call('admm2', c(list(lambda=lam), args)))
  p <- ncol(SigmaO)
  I    <- diag(p)
  args$opts <- c(args$opts, list(I=I))
  if (is.null(args$opts$tol)) args$opts$tol <- 1e-3
##  lest <- vector('list', n)
  loglik <- vector('numeric', n)
  path <- vector('list', n) ; icov <- vector('list', n) ; resid <- vector('list', n)
  for (i in n:1) {
    est <- do.call('admm2', c(list(lambda=lambda[i]), args))
    tmp <- est$S
    icov [[i]] <- tmp
    resid[[i]] <- est$L
    tmp <- Matrix::forceSymmetric(Matrix::triu(tmp,k=1))
    path [[i]] <- as(tmp, 'lsCMatrix')
    args$opts$Lambda <- est$Lambda
    args$opts$Y      <- est$Y
    if (is.null(args$opts$warm.tol)) args$opts$warm.tol <- args$opts$tol
    args$opts$tol <- args$opts$warm.tol
    z <- Matrix::rowSums(path[[i]])!=0 #1:p #
    q <- sum(!z) #p #
    R <- icov[[i]] #- resid[[i]]
    loglik[[i]] <- log(Matrix::det(R[z,z])) - sum(Matrix::diag(R[z,z] %*% SigmaO[z,z])) - (p-q)
#    args$opts$eta <- max(.5, (p-(n-i))/p)
  }
  list(icov=icov, path=path, resid=resid, lambda=lambda, loglik=loglik, data=data)
}

#' @useDynLib SpiecEasi
#' @noRd
admm2 <- function(SigmaO, lambda, beta, r, tol=1e-2, shrinkDiag=TRUE, opts) {
  n  <- nrow(SigmaO)
  defopts <- list(mu=1, eta=1, muf=1e-4, maxiter=1000, newtol=1e-4)
  if (!missing(opts)) for (o in names(opts)) defopts[[ o ]] <- opts [[ o ]]
  if (missing(beta)) beta <- 0
  if (missing(r))       r <- 0
  opts <- defopts
  if (!is.null(opts[[ 'tol' ]])) tol <- opts$tol
  over_relax_par <- 1.6
  if (!is.null(opts[[ 'over_relax_par' ]])) over_relax_par <- opts$over_relax_par
  if (!is.null(opts[[ 'over.relax.par' ]])) over_relax_par <- opts$over.relax.par

  if (is.null(opts[[ 'I' ]]))
    I <- diag(n)
  else I <- opts$I
  if (is.null(opts[[ 'Lambda' ]]))
    Lambda <- matrix(0, n, n*3)
  else Lambda <- opts$Lambda
  if (is.null(opts[[ 'Y' ]]))
    Y <- cbind(I, I, matrix(0, n, n))
  else Y <- opts$Y
  ADMM(SigmaO=SigmaO, lambda=lambda, I=I, Lambda=Lambda, Y=Y, beta=beta, r=r, shrinkDiag=shrinkDiag,
       maxiter=opts$maxiter, mu=opts$mu, eta=opts$eta, newtol=opts$newtol, muf=opts$muf,
       tol=tol, over_relax_par=over_relax_par)
}

#' robust PCA
#'
#' Form a robust PCA from clr-transformed data and \[the low rank component of\] an inverse covariance matrix
#'
#' @param X the n x p \[clr-transformed\] data
#' @param L the p x p rank-r ('residual') inverse covariance matrix from \code{spiec.easi} run argument method='slr'.
#' @param inverse flag to indicate the L is the inverse covariance matrix
#' @returns a named list with n x r matrix of scores and r x r matrix of loadings
#' @examples
#' # Create sample data
#' data(amgut1.filt)
#' data.clr <- t(clr(t(amgut1.filt), 1))
#' # Perform robust PCA
#' pca_result <- robustPCA(data.clr, diag(ncol(data.clr)))
#' @export
robustPCA <- function(X, L, inverse=TRUE) {
  Lsvd <- svd(L)
  ind <- Lsvd$d>1e-9
  if (inverse) {
    loadings <- diag(sqrt(1/Lsvd$d[ind])) %*% t(Lsvd$v[,ind])
  } else {
    loadings <- diag(sqrt(Lsvd$d[ind])) %*% t(Lsvd$v[,ind])
  }

  scores <- X %*% t(loadings)
  return(list(scores=scores, loadings=loadings))
}