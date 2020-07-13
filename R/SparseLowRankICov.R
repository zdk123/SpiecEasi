#' Sparse plus Low Rank inverse covariance
#'
#' Select an inverse covariance matrix that is a sparse plus low rank decomposition.
#'
#' @param data the n x p data matrix
#' @param npn flag to first fit nonparametric normal transform to the data
#' @param verbose flag to turn on verbose output
#' @param cor flag to use correlation matrix as the input (default: false - uses covariance)
#' @param ... arguments to override default algorithm settings (see details)
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
  args$opts$tol <- 1e-3
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
  ##  lest[[i]] <- est
    args$opts$tol <- 1
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
  defopts <- list(mu=n, eta=75/100, muf=1e-4, maxiter=100, newtol=1e-4)
  if (!missing(opts)) for (o in names(opts)) defopts[[ o ]] <- opts [[ o ]]
  if (missing(beta)) beta <- 0
  if (missing(r))       r <- 0
  opts <- defopts
  over_relax_par <- 1.6

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
       maxiter=opts$maxiter, mu=opts$mu, eta=opts$eta, newtol=opts$newtol, muf=opts$muf)
}

#' robust PCA
#'
#' Form a robust PCA from clr-transformed data and [the low rank component of] an inverse covariance matrix
#'
#' @param X the n x p [clr-transformed] data
#' @param L the p x p rank-r ('residual') inverse covariance matrix from \code{spiec.easi} run argument method='slr'.
#' @param inverse flag to indicate the L is the inverse covariance matrix
#' @returns a named list with n x r matrix of scores and r x r matrix of loadings
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
# #' @importFrom Matrix sparseMatrix forceSymmetric
# admm <- function(SigmaO, lambda, beta, r, LPD=FALSE, eigsolve=FALSE, tol=1e-5, shrinkDiag=TRUE, opts) {
#   n  <- nrow(SigmaO)
#   defopts <- list(mu=n, eta=99/100, muf=1e-4, maxiter=500, newtol=1e-4)
#   if (!missing(opts)) for (o in names(opts)) defopts[[ o ]] <- opts [[ o ]]
#   opts <- defopts
#   ABSTOL <- tol
#   RELTOL <- tol
#   over_relax_par <- 1.6
#   Ip <- sparseDiag(n)
#   Id <- diag(n)
#   J <- matrix(1/n, n, n)
#   R <- S <- Ip
#   L  <- zeros(n,n)
#   RY <- Id; SY <- Id; LY <- L
#   Y <- cbind(RY, SY, LY)
#   Lambda1 <- zeros(ncol(R)); Lambda2 <- Lambda1; Lambda3 <- Lambda1
#   mu  <- opts$mu
#   eta <- opts$eta
#   objval   <- vector('numeric', opts$maxiter)
#   r_norm   <- vector('numeric', opts$maxiter)
# #  s_norm   <- vector('numeric', opts$maxiter)
#   eps_pri  <- vector('numeric', opts$maxiter)
# #  eps_dual <- vector('numeric', opts$maxiter)
#   for (iter in 1:opts$maxiter) {
#     ## update X = (R,S,L)
#     RA <- RY + mu*Lambda1
#     SA <- SY + mu*Lambda2
#     LA <- LY + mu*Lambda3
#     if (eigsolve) {
#       tmp  <- mu*SigmaO-RA
#       tmp  <- (tmp+t(tmp))/2
#       eV   <- eigen(tmp) ; U <- eV$vectors ; D <- eV$values ; d = sparseDiag(D)
#       eigR <- (-d + sqrt(d^2+4*mu))/2
#       R    <- U%*%sparseDiag(Matrix::diag(eigR))%*%t(U)
#       R <- (R+t(R))/2
#     } else {
#       K  <- mu*SigmaO - RA
#       KI <- sqrtmNewt(as(K %*% K, 'matrix') + 4*mu*Id, Id, errTol=opts$newtol)
#       R  <- MATAVE2(KI, as.matrix(-K)) #.5*(KI-K)
#     }
#     S  <- SOFTTHRESH(as(SA, 'matrix'), lambda*mu, shrinkDiag=shrinkDiag)
# #    S  <- forceSymmetric(S)
# #    LA <- forceSymmetric(LA)
#     if (LPD) {
# #      if (iter>1) {
#         eV   <- eigen(LA) ; U <- Re(eV$vectors) ; d <- Re(eV$values)
#         if (missing(beta))
#           beta <- d[r+1]/mu #*1.01
#         eigL <- sparseDiag(pmax(d-mu*beta,0))
#         L <- U %*% eigL %*% t(U)
#         L <- as(forceSymmetric(L), 'sparseMatrix')
# #      } else L <- LA
#     } else {
#       if (iter==1) {
#         L <- LA
#       } else {
#         if (missing(beta)) {
#           tmp  <- softSVT(as(LA, 'matrix'), k=r)
# #          beta <- tmp$tau/mu
#           L    <- tmp$M
#         } else {
#             tmp <- softSVT(as(LA, 'matrix'), tau=mu*beta)
#             L   <- tmp$M
#         }
#       }
#     }
#     X  <- cbind(R, S, L)
#     RO <- over_relax_par*R + (1-over_relax_par)*RY
#     SO <- over_relax_par*S + (1-over_relax_par)*SY
#     LO <- over_relax_par*L + (1-over_relax_par)*LY
#     ## update Y = (RY,SY,LY)
#     Y_old <- Y
#     RA <- RO - mu*Lambda1
#     SA <- SO - mu*Lambda2
#     LA <- LO - mu*Lambda3
#     TA <- (RA-SA+LA)/3
#     RY <- RA - TA ; SY <- SA + TA ; LY <- LA - TA
#
#     ## update Lambda
#     Lambda1 <- (Lambda1 - (RO-RY)/mu); # Lambda1 = forceSymmetric(Lambda1)
#     Lambda2 <- (Lambda2 - (SO-SY)/mu); #Lambda2 = forceSymmetric(Lambda2)
#     Lambda3 <- (Lambda3 - (LO-LY)/mu); #Lambda3 = forceSymmetric(Lambda3)
#     Lambda  <- cbind(Lambda1, Lambda2, Lambda3)
#
#     Y <- cbind(RY, SY, LY)
#     ## diagnostics, reporting, termination checks
# #    k = iter
# #    objval[iter]  = objective(R,SigmaO,eigR,S,eigL,lambda,beta);
#
#     r_norm[iter] <- Matrix::norm(X - Y, 'F')
# #    s_norm[iter]  = max(svd(-(Y - Y_old)/mu, nu=0, nv=0)$d)
#     eps_pri[iter] <- sqrt(3*n*n)*ABSTOL +
#           RELTOL*max(Matrix::norm(X,'F'), Matrix::norm(Y,'F'))
# #    eps_dual[iter]= sqrt(3*n*n)*ABSTOL + RELTOL*Matrix::norm(Lambda,'F')
# #    print(r_norm[iter] - eps_pri[iter])
#     if (r_norm[iter] < eps_pri[iter]) # & s_norm[k] < eps_dual[k])
#       break
#   ## TODO: make this optional
#     mu <- max(mu*eta, opts$muf)
#   }
#
#   history <- list(objval=objval, r_norm=r_norm,eps_pri=eps_pri)
# # eigR = eigR, eigL = eigL,
#   list(R=R, S = S, L = L, RA=RA, LA=LA, SA=SA, mu=mu, obj = objval[iter], iter = iter, history=history)
# }
#
# zeros <- function(n, p=n, sparse=TRUE) {
#   if (!sparse)
#     matrix(0, nrow=n, ncol=p)
#   else
#     Matrix::sparseMatrix(i=NULL, j=NULL, dims=c(n,p), symmetric=TRUE)
# }
#
# sparseDiag <- function(x) {
#   if (length(x) == 1 && as.integer(x) == x) {
#     M <- zeros(x, sparse=TRUE)
#     diag(M) <- 1
#   } else {
#     M <- zeros(length(x), sparse=TRUE)
#     diag(M) <- x
#   }
#   as(M, 'dsCMatrix')
# }
#
#
# .sqrtNewton <- function(C, sqrtC0=diag(ncol(C)), errTol=5e-1) {
# # compute square root of symmetric matrix C
#     # init
#     n  <- ncol(C)
#     X <- sqrtC0
#     err <- Inf
#     while (err > errTol) {
#         X_new <- .5*(X + Matrix::solve(X, C)) #.5*sqrt(mu)*Id)))
#         err   <- Matrix::norm(X_new - X, 'F')
#         X     <- X_new
#     }
#     X
# }
#
#
# objective <- function(R,SigmaO,eigR,S,eigL,lambda,beta) {
#   sum(sum(R*SigmaO)) - sum(log(eigR)) + lambda*sum(abs(S)) + beta*sum(eigL)
# }
#
#
# hardsvdthresh <- function(M, tau, k=ncol(M)) {
#     Msvd <- svdPow(M, min(ncol(M), k+1), k+2)
# #    Msvd <- svd(M)
#     if (k < ncol(M))
#       tau <- Msvd$d[k+1]*1.01 #min(tau, Msvd$d[k+1])
#     ind <- which(Msvd$d - tau > 0)
#     tmpd <- sparseDiag(Msvd$d[ind])
#     M <- Msvd$u[,ind,drop=FALSE] %*% tmpd %*% t(Msvd$v[,ind,drop=FALSE])
#     return(list(M=M, tau=tau, d=tmpd))
# }
#
#
# svdPow <- function(A, k, q) {
#     l <- k
#     m <- nrow(A)
#     n <- ncol(A)
#     P <- matrix(rnorm(m*l), nrow=m, ncol=l)
#     QR <- qr(t(t(P)%*%A))
#     Q <- qr.Q(QR)
#     R <- qr.R(QR)
# #    return(list(Q=Q, R=R))
#     for (j in 1:q) {
#         PR <- qr(A%*%Q)
#         P  <- qr.Q(PR)
#         R  <- qr.R(PR)
#         QR <- qr(t(t(P)%*%A))
#         Q  <- qr.Q(QR)
#         R  <- qr.R(QR)
#     }
#     Asvd <- svd(A%*%Q)
#     list(u=Asvd$u, #[,1:k,drop=FALSE],
#          d=Asvd$d, #[1:k,drop=FALSE],
#          v=(Q%*%Asvd$v), P=P, Q=Q) #[,1:k,drop=FALSE])
# }
