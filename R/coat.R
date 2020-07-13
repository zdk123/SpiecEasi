#' COAT
#'
#' Compositional-adjusted thresholding by  doi.org/10.1080/01621459.2018.1442340 by Cao, Lin & Li (2018)
#' @param data a clr-transformed data or covariance matrix
#' @param lambda threshold parameter(s)
#' @param thresh "soft" or "hard" thresholding
#' @param adaptive use adative-version of the lambda as in the original COAT paper. See details.
#' @param shrinkDiag flag to exclude the covariance diagonal from the shrinkage operation
#' @param ret.icov flag to also return the inverse covariance matrix (inverse of all thresholded COAT matrices)
#' @param ... Arguments to automatically calculating the lambda path. See details.
#' @details If \code{adaptive=TRUE}, and \code{data} is a covariance matrix, the adaptive penalty is calculated by assuming the underlying data is jointly Gaussian in the infinite sample setting. The results may differ from the 'empirical' adaptive setting.
#'
#' There are a few undocumented arguments useful for computing a lambda path on the fly:
#' \itemize{
#'  \item{lambda.max: }{Maximum lambda. Default: max absolute covariance}
#'  \item{lambda.min.ratio: }{lambda.min is lambda.min.ratio*lambda.max is the smallest lambda evaluated. Defailt: 1e-3}
#'  \item{nlambda: }{Number of values of lambda between lambda.max and lambda.min. Default: 30}
#'}
#' @export
coat <- function(data, lambda, thresh="soft", adaptive=TRUE, shrinkDiag=TRUE, ret.icov=FALSE, ...) {

    if (isSymmetric(data)) {
        S <- data
        adapt <- 2
    } else {
        S <- cov(data)
        adapt <- 1
    }
    if (adaptive) {
      theta <- switch(adapt,
                '1'=.getThetaMat(data),
                '2'=sqrt(S^2 + diag(S)%*%t(diag(S))) ## assume joint gaussian model
               )
    } else
        theta <- 1

    args <- list(...)
    if (missing(lambda)) {
      if (!is.null(args[[ "lambda.max" ]])) maxlam <- args$lambda.max
      else maxlam <- pulsar::getMaxCov(abs(S)/theta)
      if (is.null(args[[ "lambda" ]])) {
        if (is.null(args[[ "lambda.min.ratio" ]])) args$lambda.min.ratio <- 1e-3
        if (is.null(args[[ "nlambda" ]])) args$nlambda <- 10
        lambda <- pulsar::getLamPath(maxlam, maxlam*args$lambda.min.ratio, args$nlambda, log=FALSE)
        args$lambda.min.ratio <- NULL ; args$nlambda <- NULL
      } else lambda <- args$lambda
    }
    args$lambda.min.ratio <- args$nlambda <- args$lambda <- args$lambda.max <- NULL



    threshfun <-
      switch(thresh,
             soft  = soft.thresh,
             hard  = hard.thresh,
             adapt = adaptive.thresh)


    if (length(lambda) > 1) {
        n <- length(lambda)
        path <- vector('list', n)
        cov  <- vector('list', n)
        if (ret.icov) icov <- vector('list', n)
        Stmp <- S
        for (i in n:1) {
            cov[[i]]  <- threshfun(Stmp, theta*lambda[i], shrinkDiag)
            if (ret.icov) icov[[i]] <- tryCatch(solve(cov[[i]]), error=function(e) MASS::ginv(cov[[i]]))
            path[[i]] <- spMat2Adj(cov[[i]])
            Stmp <- path[[i]]*S
            Matrix::diag(Stmp) <- diag(S)
        }
        est <- list(path=path, cov=cov)
    } else {
        Stmp <- threshfun(S, theta*lambda, shrinkDiag)
        if (ret.icov) icov <- tryCatch(solve(Stmp), error=function(e) MASS::ginv(Stmp))
        est <- list(
          cov  = Stmp,
          path = spMat2Adj(Stmp)
        )
    }
    est$lambda <- lambda
    if (ret.icov) est$icov <- icov
    est
}

#' @noRd
spMat2Adj <- function(A, class="lsCMatrix") {
    if (inherits(A, "sparseMatrix")) {
      return(as(Matrix::forceSymmetric(Matrix::triu(A, k=1)), class))
    } else if (inherits(A, "matrix")) {
      return(as(Matrix::forceSymmetric(Matrix::triu(A, k=1)), 'lsyMatrix'))
    } else
      stop("Class is not recognized")
}

#' @noRd
soft.thresh <- function(S, lam, shrinkDiag=FALSE) {
    if (inherits(S, "sparseMatrix")) {
        return(.sparseThresh(S, lam, shrinkDiag, thresh="soft"))
    } else if (inherits(S, "matrix")) {
        Sign <- sign(S)
        M <- as(Sign*pmax(abs(S) - lam, 0), "sparseMatrix")
        if (!shrinkDiag) Matrix::diag(M) <- diag(S)
        return(M)
    } else if (inherits(S, 'denseMatrix')) {
        return(soft.thresh(as(S, 'matrix'), lam, shrinkDiag))
    } else
        stop('Class not recognized')
}

#' @noRd
hard.thresh <- function(S, lam, shrinkDiag=FALSE) {
    if (inherits(S, "matrix")) {
        M <- Matrix::Matrix(S*(abs(S) >= lam), sparse=TRUE)
        if (!shrinkDiag) diag(M) <- diag(S)
        return(M)
    } else if (inherits(S, "sparseMatrix")) {
        return(.sparseThresh(S, lam, shrinkDiag, thresh="hard"))
    } else if (inherits(S, 'denseMatrix')) {
        return(hard.thresh(as(S, 'sparseMatrix'), lam, shrinkDiag))
    } else
        stop('Class not recognized')
}

#' @noRd
adaptive.thresh <- function(S, lam, shrinkDiag=FALSE, eta=1) {
    if (inherits(S, "matrix")) {
        p <- ncol(S)
        eta  <- pmax(eta, 1)
        M <- as(S*pmax(1 - abs(lam/S)^eta, 0), "sparseMatrix")
        if (!shrinkDiag) Matrix::diag(M) <- diag(S)
        return(M)
    } else if (inherits(S, "sparseMatrix")) {
        return(.sparseThresh(S, lam, shrinkDiag, thresh="adapt", eta))
    } else if (inherits(S, 'denseMatrix')) {
        return(adaptive.thresh(as(S, 'sparseMatrix'), lam, shrinkDiag, eta))
    } else
        stop('Class not recognized')
}

#' @noRd
.sparseThresh <- function(S, lam, shrinkDiag, thresh="soft", eta=1) {
    k <- if (shrinkDiag) 0 else 1
    p <- ncol(S)
    ilist <- Matrix::summary(Matrix::triu(S, k=k))
    ind   <- (ilist[,2]-1)*p + ilist[,1]
    if (inherits(lam, "matrix"))
        Lam <- lam[ind]
    else
        Lam <- lam
    newx <- switch(thresh,
        soft = pmax(abs(ilist[,3]) - Lam, 0),
        hard = S[ind]*(abs(ilist[,3]) >= Lam),
        adapt= S[ind]*pmax(1 - abs(Lam/ilist[,3])^eta, 0)
       )

    nzind  <- which(newx != 0)
    newx <- if (thresh=='soft') sign(ilist[nzind,3])*newx[nzind] else newx[nzind]
    M <- Matrix::sparseMatrix(ilist[nzind,1], ilist[nzind,2], x=newx, dims=dim(S))
    if (!shrinkDiag) Matrix::diag(M) <- Matrix::diag(S)
    return(Matrix::forceSymmetric(M))

}

#' @noRd
.getThetaMat <- function(X) {
  n <- ncol(X)
  theta <- matrix(0, n,n)
  for (i in 1:n) {
    for (j in i:n) {
      theta[j,i] <- theta[i,j] <- sd(X[,i]*X[,j])
    }
  }
  theta
}
