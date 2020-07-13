#' Sparse/penalized estimators of covariance matrices
#'
#' This function estimates the sparse inverse covariance matrix/matrices
#' given data (typically after clr transformation) and further arguments to
#' huge package functions.
#'
#' @param data data matrix with features/OTUs in the columns and samples in the rows. Should be transformed by clr for meaningful results, if the data is compositional
#' @param method estimation method to use as a character string. Currently either 'glasso' or 'mb' (meinshausen-buhlmann)
#' @param npn perform Nonparanormal (npn) transformation before estimation?
#' @param verbose print progress to standard out
#' @param cov.output return the covariance matrix as well.
#' @param ... further arguments to huge/estimation functions. See details.
#' @details
#' This is a wrapper function for sparse iCov estimations performed by glasso in the huge package.
#'
#' Therefore, arguments \code{...} should be named. Typically, these are for specifying a penalty parameter, \code{lambda}, or the number of penalties to use.
#' By default 10 pentalties are used, ranging logarithmically between \code{lambda.min.ratio}*MAX and MAX.
#' Max is the theoretical upper bound on \code{lambda} and us \code{max|S|}, the maximum absolute value in the data correlation matrix.
#' \code{lambda.min.ratio} is 1e-3 by default. Lower values of \code{lambda} require more memory/cpu time to compute, and sometimes huge will throw an error.
#'
#' The argument \code{nlambda} determines the number of penalties - somewhere between 10-100 is usually good, depending on how the values of empirical correlation are distributed.
#' @importFrom huge huge huge.npn
#' @export
#' @examples
#' # simulate data with 1 negative correlation
#'  set.seed(10010)
#'  Sigma <- diag(50)*2
#'  Sigma[1,2] <- Sigma[2,1] <- -.9
#'  data  <- exp(rmvnorm(50, runif(50, 0, 2), Sigma))
#'
#' # normalize
#'  data.f   <- t(apply(data, 1, norm_to_total))
#'  data.clr <- t(clr(data.f, 1))
#'
#' # estimate
#'  est.clr  <- sparseiCov(data.clr, method='glasso')
#'  est.f    <- sparseiCov(data.f, method='glasso')
#'  est.log  <- sparseiCov(log(data), method='glasso')
#'
#' # visualize results
#'  par(mfrow=c(1,3))
#'  image(as.matrix(est.log$path[[3]][1:5,1:5]))
#'  image(as.matrix(est.clr$path[[3]][1:5,1:5]))
#'  image(as.matrix(est.f$path[[3]][1:5,1:5]))
sparseiCov <- function(data, method, npn=FALSE, verbose=FALSE, cov.output = TRUE, ...) {

  if (npn) data <- huge::huge.npn(data, verbose=verbose)

  args <- list(...)
  method <- switch(method, glasso = "glasso", mb = "mb", stop("Method not supported"))

  if (is.null(args$lambda.min.ratio)) args$lambda.min.ratio <- 1e-3
  est <- do.call(huge::huge, c(args, list(x=data,
                                          method=method,
                                          verbose=verbose,
                                          cov.output = cov.output)))

  ## MB betas exported in huge>=1.3.2
  # if (method %in% c('mb')) {
  #   est <- do.call(utils::getFromNamespace('huge.mb', 'huge'),
  #                  c(args, list(x=data, verbose=verbose)))
  #   est$method <- 'mb'
  #   est$data <- data
  #   est$sym  <- ifelse(!is.null(args$sym), args$sym, 'or')
  # }
  return(est)
}


#' Neighborhood net estimates
#'
#' Select a sparse inverse covariance matrix using neighborhood selection and glmnet from various exponential models.
#' @param data n x p input (pre-transformed) data
#' @param lambda the lambda path
#' @param method ising and poisson models currently supported.
#' @param ncores number of cores for distributing the model fitting
#' @param sym symmetrize the neighborhood using the 'or' (default)/'and' rule
#' @param ... further arguments to glmnet
#' @importFrom Matrix t
neighborhood.net <- function(data, lambda, method="ising", ncores=1, sym='or', ...) {
    p <- ncol(data)
    l <- length(lambda)
    args <- list(...)
      match.fun(switch(method,
        ising   ={ nbFun <- glm.neighborhood ; args$link <- 'binomial' },
        poisson ={ nbFun <- glm.neighborhood ; args$link <- 'poisson' }
        # loglin  ={ nbFun <- llgm.neighborhood}
      ))

    ### Remove 0/1 binomials/ zero variance features
    if (method=='ising') {
      mintab <- function(x) {
          xtab <- table(x)
          length(xtab) == 1 || min(xtab) == 1
      }
      zvind  <- which(apply(data, 2, mintab))
    } else {
      zvind <- which(apply(data, 2, var)==0)
    }

    estFun <- function(i) {
      betamat      <- matrix(0, p, l)
      if (!(i %in% zvind)) {
        suppressWarnings(
          out <- do.call(nbFun,
              c(list(data[,-c(i, zvind)], data[,i,drop=FALSE], lambda), args)))
        lsub <- ncol(out)
        if (lsub<l) {
          ## extend missing lambdas value if glmnet
          ## returns only larger solution ##
          out <- cbind(out, out[,rep(lsub, l-lsub)])
        }
        betamat[-c(i, zvind),] <- out
      }
      betamat
    }

    est <- parallel::mcmapply(estFun, 1:p,
                mc.cores=ncores, SIMPLIFY='array')
    beta <- vector('list', length(lambda))
    path <- vector('list', length(lambda))
    for (i in 1:dim(est)[2]) {
        tmp       <- as(est[,i,], 'dgCMatrix')
        beta[[i]] <- tmp
        tmp       <- as(tmp, 'lgCMatrix')
        path[[i]] <- if (sym == "or") sign(tmp | t(tmp)) else sign(tmp & t(tmp))
    }
    list(beta=beta, path=path)
}


#' @importFrom glmnet glmnet
#' @noRd
glm.neighborhood <- function(X, Y, lambda, link='binomial', ...) {
    return(as.matrix(Bmat <- glmnet::glmnet(X, Y, family=link, lambda=lambda, ...)$beta))
}

# #' @useDynLib SpiecEasi LPGM_neighborhood
# llgm.neighborhood <- function(X, Y, lambda, startb=0, th=1e-6, intercept=FALSE) {
#   n = nrow(X); p = ncol(X);
#   p_new = p
#   nlams <- length(lambda)
#   X <- scale(X)
#   # NOTE: here check if intercept, change X and
#   if(intercept){
#     Xorig = X;
#     X = cbind(t(t(rep(1,n))),Xorig);
#     p_new = ncol(X);
#   }
#
#   if(length(startb) == 1 & startb == 0){startb = rep(0, p_new)}
#
#   alphasin = rep(0, nlams)
#   Bmatin = matrix(0,p,nlams);
#
#   out <- .C("LPGM_neighborhood",
#             X=as.double(t(X)), Y=as.double(Y), startb=as.double(startb),
#             lambda=as.double(lambda), n=as.integer(n), p=as.integer(p_new), nlams=as.integer(length(lambda)),
#             alphas=as.double(alphasin), Bmat=as.double(Bmatin), PACKAGE="SpiecEasi")
#
#   alphas = out$alphas
#   if(is.null(out$Bmat)) {
#     Bmat = NULL
#   } else {
#     Bmat = matrix(out$Bmat, nrow=nrow(Bmatin), byrow=TRUE)
#     Bmat <- Bmat*(abs(Bmat)>th)
#   }
#   return(Bmat)
# }

dclr <- function(x) t(clr(apply(x, 1, norm_diric),2))
dclrNPN <- function(x) huge::huge.npn(t(clr(apply(x, 1, norm_diric),2)), verbose=FALSE)
