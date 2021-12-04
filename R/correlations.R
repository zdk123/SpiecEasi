#' @export
latentcor <- function(data, types = NULL) {
  ## remove custom 'count' types ## TODO: should this be done here??
  types <- gsub("_count", "", types)
  suppressMessages(
    latentcor::latentcor(X = data, types = types, use.nearPD=FALSE, nu = 0, showplot = FALSE)$Rpointwise
  )
}

#' @keywords internal
.match.cov <- function(cov.fun, data, types) {
  switch(cov.fun,
                cor=cor(data),
                cov=cov(data),
                latentcor=latentcor(data, types=types))
}
