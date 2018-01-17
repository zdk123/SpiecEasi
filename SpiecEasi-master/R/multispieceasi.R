#' Spiec-Easi pipeline for multiple datasets.
#' All data matrices must be in a list
#' written by Laura Tipton; 1/5/18
#' @export

multi.spiec.easi <- function(datalist, method='glasso', sel.criterion='stars', verbose=FALSE, icov.select=TRUE, icov.select.params=list(), ...) {
  args <- list(...)
  if (verbose) message("Normalizing/clr transformation of each dataset with pseudocount ...")
  for (d in 1:length(datalist)){
    temppc <- datalist[[d]]+1
    temptss <- t(apply(temppc, 1, norm_to_total))
    tempclr <- t(clr(temptss,1))
    assign(paste0("clr", d), tempclr)
    #clean up intermediate matrices
    rm(temppc, temptss, tempclr)
  }
  if (verbose) message(paste("Inverse Covariance Estimation with", method, "...", sep=" "))
  est <- sparseiCov(get(paste0("crl", c(1:length(datalist)))), method=method, args)
  if (icov.select) {
    if (verbose) message(paste("Model selection with", sel.criterion, "...", sep=" "))
    sel <- do.call('icov.select', c(list(est=est, criterion=sel.criterion), icov.select.params))
    if (verbose) message("Done!")
  }
  return(sel)
}
    