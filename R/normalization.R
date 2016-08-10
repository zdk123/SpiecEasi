#####################################################################
# Different normalization schemes for microbiome counts (real or fake)
#
# @author Zachary Kurtz
# @date 10/10/2013
#####################################################################

#' add pseudocount before normalizing a vector
#' @export
norm_pseudo  <- function(x) norm_to_total(x+1)


#' Total sum normalization of a (presumably positive) data vector
#' @export
norm_to_total <- function(x) x/sum(x)


#' compute the shannon entropy from a vector (normalized internally)
#'
#' Shannon entropy is:
#'     sum [ x_i log(x_i) ]
#'      
#' @param x data vector
#' @return shannon entropy in base e
#' @export
shannon <- function(x) {
    x.f <- norm_to_total(x)
    -sum(x.f*log(x.f), na.rm=TRUE)
}

#' N_effective: Compute the exponential of the shannon entropy. linearizes shannon entropy, for a better diveristy metric (effective number of species)
#'
#' @param x data vector
#' @return N_eff in base e
#' @export
neff <- function(x) exp(shannon(x))


#' Centered log-ratio functions
#' @export
clr <- function(x, ...) {
    UseMethod('clr', x)
}

#' @method clr default
#' @export
clr.default <- function(x.f, base=exp(1), tol=.Machine$double.eps) {
    nzero <- (x.f >= tol)
    LOG <- log(ifelse(nzero, x.f, 1), base)
    ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
}

#' @method clr matrix
#' @export
clr.matrix <- function(x.f, mar=2, ...) {
    apply(x.f, mar, clr, ...)
}

#' @method clr data.frame
#' @export
clr.data.frame <- function(x.f, mar=2, ...) {
    clr(as.matrix(x.f), mar, ...)
}

#' Additive log-ratio functions
#' @export
alr <- function(x, ...) {
    UseMethod("alr", x)
}

#' @method alr default
#' @export
alr.default <- function(x.f, divcomp=1, base=exp(1), removeDivComp=TRUE,
                        tol=.Machine$double.eps) {
    zero <- (x.f >= tol)
    LOG <- log(ifelse(zero, x.f, 1), base)
    x.alr <- ifelse(zero, LOG - LOG[divcomp], 0.0)
    if (removeDivComp) x.alr[-divcomp]
    else x.alr
}


#' @method alr matrix
#' @export
alr.matrix <- function(x.f, mar=2, divcomp=1, base=exp(1), removeDivComp=TRUE,
                        tol=.Machine$double.eps) {
    if (mar == 1) x.f <- t(x.f)
    zero <- (x.f >= tol)
    LOG <- log(ifelse(zero, x.f, 1), base)
    x.alr <- ifelse(zero, LOG - LOG[,divcomp], 0.0)
    if (removeDivComp) x.alr[,-divcomp]
    else x.alr
}

#' @method alr data.frame
#' @export
alr.data.frame <- function(x.f, mar=2, ...) {
    alr(as.matrix(x.f), mar, ...)
}

#' @export
triu <- function(x) x[upper.tri(x)]
#' @export
tril <- function(x) x[lower.tri(x)]

#' @export
triu2diag <- function(x, diagval=0) {
    e <- length(x)
    n <- .5 * (sqrt(8*e + 1)+1)
    mat <- matrix(0, n, n)
    mat[upper.tri(mat)] <- x
    mat <- mat + t(mat)
    diag(mat) <- diagval
    mat
}
