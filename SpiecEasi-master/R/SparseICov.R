#' Spiec-Easi pipeline
#' Run the whole analysis, from data transformation, iCov estimation and model selection.
#' Inputs are a non-normalized OTU table and pipeline options.
#' @export
spiec.easi <- function(obj, ...) {
  UseMethod('spiec.easi', obj)
}

#' Spiec-Easi pipeline
#' @method spiec.easi phyloseq
#' @import phyloseq
#' @export
spiec.easi.phyloseq <- function(obj, ...) {
  OTU <- otu_table(obj)@.Data
  if (otu_table(obj)@taxa_are_rows) OTU <- t(OTU)
  spiec.easi.default(OTU, ...)
}


#' Spiec-Easi pipeline
#' @param data non-normalized count OTU/data table with samples on rows and features/OTUs in columns
#' @param method estimation method to use as a character string. Currently either 'glasso' or 'mb' (meinshausen-buhlmann)
#' @param sel.criterion character string specifying criterion/method for model selection accepts 'stars' [default], 'ric', 'ebic'
#' @param icov.select.params list of further arguments to icov.select
#' @param ... further arguments to sparseiCov
#' @method spiec.easi default
#' @export
spiec.easi.default <- function(data, method='glasso', sel.criterion='stars', verbose=TRUE, 
                               icov.select=TRUE, icov.select.params=list(), ...) {
  
  args <- list(...)
  if (verbose) message("Normalizing/clr transformation of data with pseudocount ...")
  data.clr <- t(clr(data+1, 1))
  if (verbose) message(paste("Inverse Covariance Estimation with", method, "...", sep=" "))
  est      <- do.call('sparseiCov', c(list(data=data.clr, method=method), args))
  
  if (icov.select) {
    if (verbose) message(paste("Model selection with", sel.criterion, "...", sep=" "))
    est <- do.call('icov.select', c(list(est=est, criterion=sel.criterion), icov.select.params))
    if (verbose) message("Done!")
  }
  return(est)
}



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
#'  Theta <- matrix(0, 50, 50)
#'  Theta[1,2] <- Theta[2,1] <- 1
#'  Sigma <- -.45*Theta
#'  diag(Sigma) <- 1
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
#' # evaluate results
#' huge::huge.roc(est.clr$path, Theta)
#' huge::huge.roc(est.log$path, Theta)
#' huge::huge.roc(est.f$path,   Theta)
#'
sparseiCov <- function(data, method, npn=FALSE, verbose=FALSE, cov.output = TRUE, ...) {
  
  if (npn) data <- huge::huge.npn(data, verbose=verbose)

  args <- list(...)

  method <- switch(method, glasso = "glasso", mb = "mb", stop("Method not supported"))

  if (is.null(args$lambda.min.ratio)) args$lambda.min.ratio <- 1e-3

  if (method %in% c("glasso")) {
    do.call(huge::huge, c(args, list(x=data, method=method, verbose=verbose, 
                                     cov.output = cov.output)))

  } else if (method %in% c('mb')) {
    est <- do.call(huge::huge.mb, c(args, list(x=data, verbose=verbose)))
    est$method <- 'mb'
    est$data <- data
    est$sym  <- ifelse(!is.null(args$sym), args$sym, 'or')
    return(est)
  }
}


#' Model selection for picking the right \code{lambda} penalty.
#' This is identical to huge::huge.stars except that the subsampling loop is replaced with an mclapply function to add parallelization capabilities.
#' 
#' @param est an estimate/model as produced by the sparseiCov function
#' @param criterion character string specifying criterion/method for model selection accepts 'stars' [default], 'ric', 'ebic'
#' @param stars.thresh variability threshold for stars selection
#' @param ebic.gamma tuning parameter for ebic
#' @param stars.subsample.ratio The default value 'is 10*sqrt(n)/n' when 'n>144' and '0.8' when 'n<=144', where 'n' is the sample size.
#' @param rep.num number of subsamplings when \code{criterion} = stars.
#' @param ncores number of cores to use. Need multiple processers if \code{ncores > 1}
#' @param normfun normalize internally if data should be renormalized
#' @importFrom parallel mclapply
#' @export
icov.select <- function(est, criterion = 'stars', stars.thresh = 0.05, ebic.gamma = 0.5, 
                        stars.subsample.ratio = NULL, rep.num = 20, ncores=1, normfun=function(x) x, verbose=FALSE) {
  gcinfo(FALSE)
  if (est$cov.input) {
    message("Model selection is not available when using the covariance matrix as input.")
    class(est) = "select"
    return(est)
  }
  if (!est$cov.input) {
    if (est$method == "mb" && is.null(criterion)) 
      criterion = "stars"
    if (est$method == "ct" && is.null(criterion)) 
      criterion = "ebic"
    n = nrow(est$data)
    d = ncol(est$data)
    nlambda = length(est$lambda)
    if (criterion == "ric") {
      if (verbose) {
        message("Conducting rotation information criterion (ric) selection....")
#        flush.console()
      }
      if (n > rep.num) {
        nr = rep.num
        r = sample(n, rep.num)
      }
      if (n <= rep.num) {
        nr = n
        r = 1:n
      }
      out = .C("RIC", X = as.double(est$data), dd = as.integer(d), 
               nn = as.integer(n), r = as.integer(r), nr = as.integer(nr), 
               lambda_opt = as.double(0), PACKAGE = "huge")
      est$opt.lambda = out$lambda_opt/n
      rm(out)
      gc()
      if (verbose) {
        message("done\n")
#        flush.console()
      }
      if (verbose) {
        message("Computing the optimal graph....")
#        flush.console()
      }
      if (est$opt.lambda > max(cor(est$data))) 
        est$refit = Matrix(0, d, d)
      else {
        if (est$method == "mb") 
          est$refit = huge::huge.mb(est$data, lambda = est$opt.lambda, 
                              sym = est$sym, idx.mat = est$idx.mat, verbose = FALSE)$path[[1]]
        if (est$method == "glasso") {
          if (!is.null(est$cov)) {
            tmp = huge::huge.glasso(est$data, lambda = est$opt.lambda, 
                              scr = est$scr, cov.output = TRUE, verbose = FALSE)
            est$opt.cov = tmp$cov[[1]]
          }
          if (is.null(est$cov)) 
            tmp = huge::huge.glasso(est$data, lambda = est$opt.lambda, 
                              verbose = FALSE)
          est$refit = tmp$path[[1]]
          est$opt.icov = tmp$icov[[1]]
          rm(tmp)
          gc()
        }
        if (est$method == "ct") 
          est$refit = huge::huge.ct(est$data, lambda = est$opt.lambda, 
                              verbose = FALSE)$path[[1]]
      }
      est$opt.sparsity = sum(est$refit)/d/(d - 1)
      if (verbose) {
        cat("done\n")
#        flush.console()
      }
    }
    if (criterion == "ebic" && est$method == "glasso") {
      if (verbose) {
        cat("Conducting extended Bayesian information criterion (ebic) selection....")
#        flush.console()
      }
      est$ebic.score = -n * est$loglik + log(n) * est$df + 4 * ebic.gamma * log(d) * est$df
      est$opt.index = which.min(est$ebic.score)
      est$refit = est$path[[est$opt.index]]
      est$opt.icov = est$icov[[est$opt.index]]
      if (est$cov.output) 
        est$opt.cov = est$cov[[est$opt.index]]
      est$opt.lambda = est$lambda[est$opt.index]
      est$opt.sparsity = est$sparsity[est$opt.index]
      if (verbose) {
        message("done\n")
#        flush.console()
      }
    }
    if (criterion == "stars") {
      if (is.null(stars.subsample.ratio)) {
        if (n > 144) 
          stars.subsample.ratio = 10 * sqrt(n)/n
        if (n <= 144) 
          stars.subsample.ratio = 0.8
      }
      
      #            for (i in 1:nlambda) merge[[i]] <- Matrix(0, d, d)
      
      if (verbose) {
        mes = "Conducting Subsampling....."
        message(mes, appendLF = FALSE)
#        cat("\n")
#        flush.console()
      }
      #    for (i in 1:rep.num) {
      premerge <- parallel::mclapply(1:rep.num, function(i) {
        #                if (verbose) {
        #                  mes <- paste(c("Conducting Subsampling....in progress:", 
        #                    floor(100 * i/rep.num), "%"), collapse = "")
        #                  cat(mes, "\r")
        #                  flush.console()
        #                }
        #                merge <- replicate(nlambda, Matrix(0, d,d))
        ind.sample = sample(c(1:n), floor(n * stars.subsample.ratio), 
                            replace = FALSE)
        if (est$method == "mb") 
          tmp = huge::huge.mb(normfun(est$data[ind.sample, ]), lambda = est$lambda, 
                        scr = est$scr, idx.mat = est$idx.mat, sym = est$sym, 
                        verbose = FALSE)$path
        if (est$method == "ct") 
          tmp = huge::huge.ct(normfun(est$data[ind.sample, ]), lambda = est$lambda, 
                        verbose = FALSE)$path
        if (est$method == "glasso") 
          tmp = huge::huge.glasso(normfun(est$data[ind.sample, ]), lambda = est$lambda, 
                            scr = est$scr, verbose = FALSE)$path
        #                for (j in 1:nlambda) merge[[j]] <- merge[[j]] + tmp[[j]]

        rm(ind.sample)
        gc()
        return(tmp)
      }, mc.cores=ncores)
      #  }
      # merge <- lapply(merge, as.matrix)
#      merge <- lapply(merge, simplify2array)
#      est$merge <- lapply(1:dim(merge)[3], function(i) merge[,,i]/rep.num)

      merge <- Reduce(function(l1, l2) lapply(1:length(l1),
                    function(i) l1[[i]] + l2[[i]]), premerge, accumulate=FALSE)

      if (verbose) {
        message("done")
#        cat("\n")
#        flush.console()
      }
      est$variability = rep(0, nlambda)
      est$merge <- vector('list', nlambda)
      for (i in 1:nlambda) {
        est$merge[[i]]  <- merge[[i]]/rep.num
        est$variability[i] <- 4 * sum(est$merge[[i]] * (1 - est$merge[[i]]))/(d * (d - 1))
      }
      est$opt.index = max(which.max(est$variability >= 
                                      stars.thresh)[1] - 1, 1)
      est$refit = est$path[[est$opt.index]]
      est$opt.lambda = est$lambda[est$opt.index]
      est$opt.sparsity = est$sparsity[est$opt.index]
      if (est$method == "glasso") {
        est$opt.icov = est$icov[[est$opt.index]]
        if (!is.null(est$cov)) 
          est$opt.cov = est$cov[[est$opt.index]]
      }
    }
    est$criterion = criterion
    class(est) = "select"
    return(est)
  }
}


#' @keywords internal
dclr <- function(x) t(clr(apply(x, 1, norm_diric),2))
#' @keywords internal
dclrNPN <- function(x) huge::huge.npn(t(clr(apply(x, 1, norm_diric),2)), verbose=FALSE)
