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

  method <- switch(method, glasso = "glasso", mb = "mb",
                   stop("Method not supported"))

  if (is.null(args$lambda.min.ratio))
    args$lambda.min.ratio <- 1e-3

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
