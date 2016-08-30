coat <- function(data, lambda, thresh="soft", adaptive=TRUE, shrinkDiag=FALSE, ...) {

    if (isSymmetric(data)) {
        S <- data
        adaptive <- FALSE
    } else
        S <- cov(data)

    if (adaptive) {
        theta <- getThetaMat(scale(data))
    } else
        theta <- 1

    threshfun <-
      switch(thresh,
             soft  = soft.thresh,
             hard  = hard.thresh,
             adapt = adaptive.thresh)


    if (length(lambda) > 1) {
        n <- length(lambda)
        path <- vector('list', n)
        cov  <- vector('list', n)
        Stmp <- S
        for (i in n:1) {
            cov[[i]]  <- threshfun(Stmp, theta*lambda[i], shrinkDiag, ...)
            path[[i]] <- spMat2Adj(cov[[i]])
            Stmp <- path[[i]]*S
            Matrix::diag(Stmp) <- diag(S)
        }
        est <- list(path=path, cov=cov)
    } else {
        Stmp <- threshfun(S, theta*lambda, shrinkDiag, ...)
        est <- list(
          cov  = Stmp,
          path = spMat2Adj(Stmp)
        )
    }
    est
}


spMat2Adj <- function(A, class="lsCMatrix") {
    if (inherits(A, "sparseMatrix")) {
      return(as(forceSymmetric(Matrix::triu(A, k=1)), class))
    } else if (inherits(A, "matrix")) {
      return(as(forceSymmetric(Matrix::triu(A, k=1)), 'lsyMatrix'))
    } else
      stop("Class is not recognized")
}


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

hard.thresh <- function(S, lam, shrinkDiag=FALSE) {
    if (inherits(S, "matrix")) {
        M <- Matrix(S*(abs(S) >= lam), sparse=TRUE)
        if (!shrinkDiag) diag(M) <- diag(S)
        return(M)
    } else if (inherits(S, "sparseMatrix")) {
        return(.sparseThresh(S, lam, shrinkDiag, thresh="hard"))
    } else if (inherits(S, 'denseMatrix')) {
        return(hard.thresh(as(S, 'sparseMatrix'), lam, shrinkDiag))
    } else 
        stop('Class not recognized')
}

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
    M <- sparseMatrix(ilist[nzind,1], ilist[nzind,2], x=newx, dims=dim(S))
    if (!shrinkDiag) Matrix::diag(M) <- Matrix::diag(S)
    return(forceSymmetric(M))

}

sdprod <- function(x, y) sd(x*y)

getThetaMat <- function(X) {
    kernlab::kernelMatrix(sdprod, t(X))@.Data
}
