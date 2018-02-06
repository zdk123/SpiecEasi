#' SPIEC-EASI pipeline
#'
#' Run the whole SPIEC-EASI pipeline, from data transformation, sparse inverse covariance estimation and model selection.
#' Inputs are a non-normalized OTU table and pipeline options.
#' @export
spiec.easi <- function(obj, ...) {
  UseMethod('spiec.easi', obj)
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

#' @param obj phyloseq object or just the otu_table
#' @method spiec.easi phyloseq
#' @rdname spiec.easi
#' @export
spiec.easi.phyloseq <- function(obj, ...) {
  if (!require('phyloseq')) {
    stop('\'Phyloseq\' package is not installed. See doi.org/doi:10.18129/B9.bioc.phyloseq')
  }
  spiec.easi.default(.phy2mat(obj), ...)
}

#' @method spiec.easi otu_table
#' @rdname spiec.easi
#' @export
spiec.easi.otu_table <- function(obj, ...) {
  if (!require('phyloseq')) {
    stop('\'Phyloseq\' package is not installed. See doi.org/doi:10.18129/B9.bioc.phyloseq')
  }
  spiec.easi.default(.phy2mat(obj), ...)

}


#' @param data non-normalized count OTU/data table with samples on rows and features/OTUs in columns
#' @param method estimation method to use as a character string. Currently either 'glasso' or 'mb' (meinshausen-buhlmann's neighborhood selection)
#' @param sel.criterion character string specifying criterion/method for model selection. Accepts 'stars' [default], 'ric', 'ebic' (not recommended for high dimensional data)
#' @param verbose flag to show progress messages
#' @param icov.select.params list of further arguments to \code{\link{icov.select}}
#' @param ... further arguments to \code{\link{sparseiCov}} / \code{huge}
#' @method spiec.easi default
#' @rdname spiec.easi
#' @seealso multi.spiec.easi
#' @export
spiec.easi.default <- function(data, method='glasso', sel.criterion='stars', verbose=TRUE,
                               icov.select=TRUE, icov.select.params=list(), ...) {

  args <- list(...)
  data.clr <- .spiec.easi.norm(data, verbose)

  if (verbose)
    message(paste("Inverse Covariance Estimation with", method, "...", sep=" "))
  est <- do.call('sparseiCov', c(list(data=data.clr, method=method), args))

  if (icov.select) {
    if (verbose) message(paste("Model selection with", sel.criterion, "...", sep=" "))
    est <- do.call('icov.select', c(list(est=est, criterion=sel.criterion), icov.select.params))
    if (verbose) message("Done!")
  }
  return(est)
}


#' @keywords internal
.spiec.easi.norm <- function(data, verbose) {
# internal function to normalize a data matrix
  if (verbose)
    message("Normalizing/clr transformation of data with pseudocount ...")

  if (inherits(data, 'matrix')) {
    ## standard data pipeline
    return(t(clr(data+1, 1)))
  } else if (inherits(data, 'list')) {
    ## multi domain spiec.easi, data must be list of numeric matrices
    return(do.call('cbind', lapply(data, .spiec.easi.norm, verbose=FALSE)))
  } else {
    stop('input data must be a numeric matrix')
  }
}

#' multi domain SPIEC-EASI
#'
#' A SPIEC-EASI pipeline for inferring a sparse inverse covariance matrix within and between multiple compositional datasets, under joint sparsity penalty.
#'
#' Can also run \code{spiec.easi} on a list and S3 will dispatch the proper function.
#' @param datalist list of non-normalized count OTU/data tables (stored in a matrix, data.frame or phyloseq/otu_table) with samples on rows and features/OTUs in columns
#' @param method estimation method to use as a character string. Currently either 'glasso' or 'mb' (meinshausen-buhlmann's neighborhood selection)
#' @param sel.criterion character string specifying criterion/method for model selection. Accepts 'stars' [default], 'ric', 'ebic' (not recommended for high dimensional data)
#' @param verbose flag to show progress messages
#' @param icov.select.params list of further arguments to \code{\link{icov.select}}
#' @param ... further arguments to \code{\link{sparseiCov}} / \code{huge}
#' @seealso spiec.easi
#' @export
multi.spiec.easi <- function(datalist, method='glasso', sel.criterion='stars', verbose=TRUE, icov.select=TRUE, icov.select.params=list(), ...) {
## functional wrapper for spiec.easi.list
  spiec.easi.list(datalist, method=method, sel.criterion=sel.criterion,
      verbose=verbose, icov.select=icov.select,
      icov.select.params=icov.select.params, ...)
}

#' @method spiec.easi list
#' @rdname multi.spiec.easi
#' @export
spiec.easi.list <- function(obj, ...) {
  classes <- sapply(obj, class)
  if (length(unique(classes)) != 1)
    warning('input list contains data of mixed classes.')

  ## convert phyloseq objects to matrices
  if (any('phyloseq' %in% classes) || any('otu_table' %in% classes))
    obj <- lapply(obj, .phy2mat)

  ## Finally, check the number of rows (samples) are equal
  ## and sample names are identical (sample names can be NULL)
  ssizes <- lapply(obj, nrow)
  snames <- lapply(obj, row.names)
  list.equal <- function(li) sum(duplicated(li)) == length(li)-1

  if (!list.equal(snames) || !list.equal(ssizes))
    stop("Do not run multi.spiec.easi with unidentical sample scheme")

  spiec.easi.default(obj, ...)
}
