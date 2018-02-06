#' Spiec-Easi pipeline
#'
#' Run the whole analysis, from data transformation, iCov estimation and model selection.
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

#' Spiec-Easi pipeline
#' @method spiec.easi phyloseq
#' @export
spiec.easi.phyloseq <- function(obj, ...) {
  if (!require('phyloseq')) {
    stop('\'phyloseq\' package is not installed')
  }
  spiec.easi.default(.phy2mat(obj), ...)
}

#' Spiec-Easi pipeline
#' @method spiec.easi list
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


#' internal function to normalize a data matrix
#' @keywords internal
.spiec.easi.norm <- function(data, verbose) {
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


#' Spiec-Easi pipeline for multiple compositional datasets.
#' All data matrices must be in a list
#' written by Laura Tipton; 1/5/18
#' @export
multi.spiec.easi <- function(datalist, method='glasso', sel.criterion='stars', verbose=TRUE, icov.select=TRUE, icov.select.params=list(), ...) {
## functional wrapper for spiec.easi.list
  spiec.easi.list(datalist, method=method, sel.criterion=sel.criterion,
      verbose=verbose, icov.select=icov.select,
      icov.select.params=icov.select.params, ...)
}
