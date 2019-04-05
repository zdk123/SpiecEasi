#' SPIEC-EASI pipeline
#'
#' Run the whole SPIEC-EASI pipeline, from data transformation, sparse inverse covariance estimation and model selection.
#' Inputs are a non-normalized OTU table and pipeline options.
#' @export
#' @importFrom pulsar pulsar batch.pulsar
spiec.easi <- function(data, ...) {
  UseMethod('spiec.easi', data)
}

#' @keywords internal
.phy2mat <- function(OTU) {
  if (inherits(OTU, 'phyloseq'))
    OTU <- OTU@otu_table
  if (inherits(OTU, 'otu_table')) {
    if (OTU@taxa_are_rows) OTU <- t(OTU@.Data)
    else OTU <- OTU@.Data
  }
  return(OTU)
}

#' @keywords internal
.data.checks <- function(data) {
  ## data checks ##
  if (inherits(data, 'list')) {
    sink <- lapply(data, .data.checks)
    return(NULL)
  }
  if (isTRUE(all.equal(rowSums(data), rep(1L, nrow(data))))) {
    warning('Data is normalized, but raw counts are expected')
  }

  if (any(data<0)) {
    warning('Negative values detected, but raw counts are expected')
  }
  return(NULL)
}

#' @method spiec.easi phyloseq
#' @rdname spiec.easi
#' @export
spiec.easi.phyloseq <- function(data, ...) {
  if (!requireNamespace('phyloseq', quietly=TRUE)) {
    stop('\'Phyloseq\' package is not installed. See doi.org/doi:10.18129/B9.bioc.phyloseq')
  }
  spiec.easi(.phy2mat(data), ...)
}

#' @method spiec.easi otu_table
#' @rdname spiec.easi
#' @export
spiec.easi.otu_table <- function(data, ...) {
  if (!requireNamespace('phyloseq', quietly=TRUE)) {
    stop('\'Phyloseq\' package is not installed. See doi.org/doi:10.18129/B9.bioc.phyloseq')
  }
  spiec.easi(.phy2mat(data), ...)

}



#' @keywords internal
.spiec.easi.norm <- function(data) {
# internal function to normalize a data matrix
  if (inherits(data, 'matrix')) {
    ## standard data pipeline
    return(t(clr(data+1, 1)))
  } else if (inherits(data, 'list')) {
    ## multi domain spiec.easi, data must be list of numeric matrices
    return(do.call('cbind', lapply(data, .spiec.easi.norm)))
  } else {
    stop('input data must be a numeric matrix')
  }
}


#' @name pulsar.params
#' @title pulsar params
#' @description The values to the \code{pulsar.params}/\code{icov.select.params} argument in the \code{\link{spiec.easi}} function must be a list with values from pulsar and/or batch.pulsar. See the pulsar docs for detailed instructions.
#'
#' List of arguments, data type, default. Description
#' \itemize{
#'    \item thresh, numeric, 0.05. Threshold for StARS criterion.
#'    \item subsample.ratio, numeric, 0.8. Subsample size for StARS.
#'    \item rep.num, numeric, 20. Number of subsamples for StARS.
#'    \item seed, numeric, NULL. Set the random seed for subsample set.
#'    \item ncores, numeric, 1. Number of cores for parallel.
#' }
#' With \code{pulsar.select='batch'}, additional arguments:
#' \itemize{
#'    \item wkdir, dir path, current directory. Working directory for process running jobs.
#'    \item regdir, dir path, temp directory. Directory for storing the registry files.
#'    \item init, string, 'init'. String for differentiating the init registry for batch mode pulsar.
#'    \item conffile, string / file path, ''. Path to config file or string that identifies a default config file.
#'    \item job.res, list, empty list. Named list to specify job resources for an hpc.
#'    \item cleanup, boolean, FALSE. Remove registry files.
#'}
#' @seealso \code{\link[pulsar]{pulsar}} \code{\link[pulsar]{batch.pulsar}} \code{\link{spiec.easi}}
NULL

#' @keywords internal
.check_pulsar_params <- function(fun, args=list()) {
  if (!inherits(args, 'list') || (length(args) >0 && is.null(names(args))) || any('' %in% names(args))) {
    stop('pulsar.params must be a named list')
  }
  if (length(args)==0) return(TRUE)
  fun   <- match.fun(fun)
  forms <- formals(fun)

  ## disallowed arguments provided by spiec.easi
  nargs    <- c("data", "fun", "fargs", "criterion")
  extrargs <- intersect(names(args), nargs)
  if (length(extrargs)>0) {
    argstr <- paste0(paste0("'", extrargs, "'"), collapse=", ")
    sstr <- ifelse(length(extrargs)==1, '', 's')
    stop(sprintf("Disallowed argument%s to \'pulsar.params\': %s", sstr, argstr))
  }
  ## check provided args
  allforms <- setdiff(names(forms), nargs)
  extrargs <- setdiff(names(args), allforms)
  if (length(extrargs)>0) {
    argstr <- paste0(paste0("'", extrargs, "'"), collapse=", ")
    sstr <- ifelse(length(extrargs)==1, '', 's')
    stop(sprintf("Unrecognized argument%s to \'pulsar.params\': %s", sstr, argstr))
  }
  return(TRUE)
}


#' @param data For a matrix, non-normalized count OTU/data table with samples on rows and features/OTUs in columns. Can also by phyloseq or otu_table object.
#' @param method estimation method to use as a character string. Currently either 'glasso' or 'mb' (meinshausen-buhlmann's neighborhood selection)
#' @param sel.criterion character string specifying criterion/method for model selection. Accepts 'stars' [default], 'bstars' (Bounded StARS)
#' @param verbose flag to show progress messages
#' @param pulsar.select flag to perform model selection. Choices are TRUE/FALSE/'batch'
#' @param pulsar.params list of further arguments to \code{\link{pulsar}} or \code{\link{batch.pulsar}}. See the documentation for \code{\link{pulsar.params}}.
#' @param icov.select deprecated.
#' @param icov.select.params deprecated.
#' @param ... further arguments to \code{\link{sparseiCov}} / \code{huge}
#' @method spiec.easi default
#' @rdname spiec.easi
#' @seealso multi.spiec.easi
#' @export
spiec.easi.default <- function(data, method='glasso', sel.criterion='stars',
                        verbose=TRUE, pulsar.select=TRUE, pulsar.params=list(),
                        icov.select=pulsar.select,
                        icov.select.params=pulsar.params, ...) {

  args <- list(...)

  .data.checks(data)

  if (verbose) msg <- .makeMessage("Applying data transformations...")
  else msg <- .makeMessage('')

  switch(method,
        glasso = {
                    message(msg, appendLF=verbose)
                    estFun <- "sparseiCov"
                    args$method <- method
                    X <- .spiec.easi.norm(data)
                    if (is.null(args[['lambda.max']]))
                      args$lambda.max <- pulsar::getMaxCov(cor(X))
                 },

        mb     = {
                    message(msg, appendLF=verbose)
                    estFun <- "sparseiCov"
                    args$method <- method
                    X <- .spiec.easi.norm(data)
                    if (is.null(args[['lambda.max']]))
                      args$lambda.max <- pulsar::getMaxCov(cor(X))
                  },

    )

  if (is.null(args[[ "lambda" ]])) {
    if (is.null(args[[ "lambda.min.ratio" ]])) args$lambda.min.ratio <- 1e-3
    if (is.null(args[[ "nlambda" ]])) args$nlambda <- 20
    args$lambda <- pulsar::getLamPath(args$lambda.max, args$lambda.max*args$lambda.min.ratio,
                              args$nlambda, log=TRUE)
    args$lambda.min.ratio <- args$nlambda <- args$lambda.max <- NULL
  }

  ocall <- match.call(expand.dots=FALSE)
  ## if pulsar options are not specified, check for deprecated icov.select options are
  if (is.null(ocall[["pulsar.select"]]) && is.null(ocall[["pulsar.params"]])) {
    pulsar.select <- icov.select
    pulsar.params <- icov.select.params
  }

  if (!is.null(pulsar.params[[ "data" ]]))
    stop("supply data directly to spiec.easi, not pulsar.params")
  if (!is.null(pulsar.params[[ "criterion" ]]))
    stop("supply sel.criterion directly to spiec.easi, not pulsar.params")


  if (pulsar.select=="batch") {
    fun <- "batch.pulsar"
    call <- quote(batch.pulsar(data=X, fun=match.fun(estFun), fargs=args))
    pulsar.select <- TRUE
  } else {
    fun <- "pulsar"
    call <- quote(pulsar(data=X, fun=match.fun(estFun), fargs=args))
  }

  if (pulsar.select) {
    ## process pulsar.params defaults
    flag <- .check_pulsar_params(fun, pulsar.params)

    pulsar.params$criterion <-
      switch(sel.criterion,
             stars = "stars",
            bstars = "stars",
  #          gstars = c("stars", "gcd"), #TODO: process gstars option
            stop("Unknown selection criterion"))

    if (sel.criterion %in% c("bstars", "gstars"))
      pulsar.params$lb.stars <- pulsar.params$ub.stars <- TRUE
    if (is.null(pulsar.params[[ "thresh" ]])) pulsar.params$thresh <- 0.05

    obj <- list(call=call)
    class(obj) <- 'pulsar'
    call <- do.call('update',
              c(pulsar.params, list(object=obj, evaluate=FALSE)))

    if (verbose)
      message(sprintf("Selecting model with %s using ", fun), sel.criterion, "...")
    est <- eval(call, environment())
    if (sel.criterion == "gstars")
      opt.index <- pulsar::opt.index(est) <- pulsar::get.opt.index(est, 'gcd')
    else
      opt.index <- pulsar::opt.index(est, 'stars')
  } else
    est <- structure(list(call=call, envir=environment()), class='pulsar')

  if (verbose) message("Fitting final estimate with ", method, "...")
  suppressWarnings(
  fit <- pulsar::refit(est)
  )
  if (pulsar.select) {
    fit$select <- est
#    list(
#     merge     = est$stars$merge,
#     summary   = est$stars$summary,
#     opt.index = opt.index,
#     lb.index  = est$stars$lb.index,
#     ub.index  = est$stars$ub.index
#    )
#     fit$refit <- if (sel.criterion=="gstars") fit$refit$gcd else fit$refit$stars
#     attr(fit$refit, 'names') <- sel.criterion
  }
  fit$lambda <- args$lambda
  fit$fun    <- call(estFun)[[1]]
  if (verbose) message('done')

  return(fit)
}


#' multi domain SPIEC-EASI
#'
#' A SPIEC-EASI pipeline for inferring a sparse inverse covariance matrix within and between multiple compositional datasets, under joint sparsity penalty.
#'
#' Can also run \code{spiec.easi} on a list and S3 will dispatch the proper function.
#' @param datalist list of non-normalized count OTU/data tables (stored in a matrix, data.frame or phyloseq/otu_table) with samples on rows and features/OTUs in columns
#' @param method estimation method to use as a character string. Currently either 'glasso' or 'mb' (meinshausen-buhlmann's neighborhood selection)
#' @param sel.criterion character string specifying criterion/method for model selection. Accepts 'stars' and 'bstars' [default]
#' @param verbose flag to show progress messages
#' @param pulsar.select flag to perform model selection. Choices are TRUE/FALSE/'batch'
#' @param pulsar.params list of further arguments to \code{\link{pulsar}} or \code{\link{batch.pulsar}}. See the documentation for \code{\link{pulsar.params}}.
#' @param ... further arguments to \code{\link{sparseiCov}} / \code{huge}
#' @seealso spiec.easi
#' @export
multi.spiec.easi <- function(datalist, method='glasso', sel.criterion='stars',
                        verbose=TRUE, pulsar.select=TRUE, pulsar.params=list(),
                        ...) {
## functional wrapper for spiec.easi.list
  spiec.easi.list(datalist, method=method, sel.criterion=sel.criterion,
                  verbose=verbose, pulsar.select=pulsar.select,
                  pulsar.params=pulsar.params, ...)
}

#' @method spiec.easi list
#' @param data non-normalized count OTU/data table with samples on rows and features/OTUs in columns. Can also be list of phyloseq objects.
#' @rdname multi.spiec.easi
#' @export
spiec.easi.list <- function(data, ...) {
  classes <- sapply(data, class)
  if (length(unique(classes)) != 1)
    warning('input list contains data of mixed classes.')

  ## convert phyloseq objects to matrices
  if (any('phyloseq' %in% classes) || any('otu_table' %in% classes))
    data <- lapply(data, .phy2mat)

  ## Finally, check the number of rows (samples) are equal
  ## and sample names are identical (sample names can be NULL)
  ssizes <- lapply(data, nrow)
  snames <- lapply(data, row.names)
  list.equal <- function(li) sum(duplicated(li)) == length(li)-1

  if (!list.equal(snames) || !list.equal(ssizes))
    stop("Do not run multi.spiec.easi with unidentical sample scheme")
  spiec.easi.default(data, ...)
}
