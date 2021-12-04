#' @export
is.count <- function(x) {
  stopifnot(is.vector(x))
  if (!is.numeric(x) && !is.integer(x)) return(FALSE)
  if (isFALSE(all.equal(x, as.integer(x)))) return(FALSE)
  if (any(x<0)) return(FALSE)
  return(TRUE)
}

#' @export
is.normalized <- function(X) {
  # assuming n x p matrix
  stopifnot(is.matrix(X) | is.data.frame(X))
  isTRUE(all.equal(rowSums(X), rep(1L, nrow(X))))
}


#' @export
get_types <- function(X, ...) {
  types <- latentcor::get_types(X)
  ## determine count and compositional types as a subset of truncated or continuous
  tc_ind <- (types=="tru") | (types=="con")
  is_count <- apply(X[,tc_ind], 2, is.count)
  types[tc_ind][is_count] <- paste(types[tc_ind][is_count], "count", sep="_")
  types
}
