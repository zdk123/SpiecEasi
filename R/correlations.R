#' @export
latentcor <- function(data, types = NULL) {
  suppressMessages(
    latentcor::latentcor(X = data, types = types, use.nearPD=FALSE, nu = 0, showplot = FALSE)$Rpointwise
  )
}
