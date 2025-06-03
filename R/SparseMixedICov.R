# sparseMixediCov <- function(data, types, method, verbose=FALSE, cov.output = TRUE, ...) {

#   # data is a list of n-row matrices of length r
#     # in which case types is a vector of length r
#   # TODO: accept preconcatenated data types?
#   stopifnot(inherits(data, 'list'))
#   stopifnot(sapply(data, inherits, 'matrix'))
#   stopifnot(length(data)==length(types))
#   if (!is.null(types)) {
#     stopifnot(types %in% c("con", "bin", "tru", "ter"))
#   }
#   ## Let latentcor infer types
#   # else {
#   #   stop("specify data types")
#   # }

# ##  estR <- mixedCCA:::.estimateR_mixed_multi(data, types, use.nearPD=FALSE, nu=0)$R
#   ## TODO: check if use of use.nearPD in this function is a problem (should be optional!?)
#   estR <- latentcor::latentcor(X = data, types = types, nu = 0)$R

#   SpiecEasi::sparseiCov(estR, method, npn, verbose, cov.output, ...)
# }
