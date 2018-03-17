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
      }
      if (verbose) {
        message("Computing the optimal graph....")
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
        message("done")
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
