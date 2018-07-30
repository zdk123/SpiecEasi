#####################################################################
# Different normalization schemes for microbiome counts (real or fake)
#
# @author Zachary Kurtz
# @date 10/10/2013
#####################################################################

#' Normalize w/ Pseudocount
#'
#' add pseudocount before normalizing a count vector
#'
#' @param x count data vector
#' @export
norm_pseudo  <- function(x) norm_to_total(x+1)

#' Total Sum Normalize
#'
#' Normalize a count vector by the total sum of that vector
#' @param x count data vector
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
#' @param x.f input data
#' @param ... pass through arguments
#' @export
clr <- function(x.f, ...) {
    UseMethod('clr')
}

#' @method clr default
#' @param base base for log transformation
#' @param tol tolerance for a numerical zero
#' @rdname clr
#' @export
clr.default <- function(x.f, base=exp(1), tol=.Machine$double.eps, ...) {
    nzero <- (x.f >= tol)
    LOG <- log(ifelse(nzero, x.f, 1), base)
    ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
}

#' @method clr matrix
#' @param mar margin to apply the transformation (rows: 1 or cols: 2)
#' @rdname clr
#' @export
clr.matrix <- function(x.f, mar=2, ...) {
    apply(x.f, mar, clr, ...)
}

#' @method clr data.frame
#' @rdname clr
#' @export
clr.data.frame <- function(x.f, mar=2, ...) {
    clr(as.matrix(x.f), mar, ...)
}

#' Additive log-ratio functions
#' @param x.f input data
#' @param ... pass through arguments
#' @export
alr <- function(x.f, ...) {
    UseMethod("alr")
}

#' @method alr default
#' @param divcomp the index of the component to use as the divisor
#' @param base base for log transformation
#' @param removeDivComp remove the divisor component from the alr result
#' @param tol tolerance for a numerical zero
#' @rdname alr
#' @export
alr.default <- function(x.f, divcomp=1, base=exp(1), removeDivComp=TRUE,
                        tol=.Machine$double.eps, ...) {
    zero <- (x.f >= tol)
    LOG <- log(ifelse(zero, x.f, 1), base)
    x.alr <- ifelse(zero, LOG - LOG[divcomp], 0.0)
    if (removeDivComp) x.alr[-divcomp]
    else x.alr
}


#' @method alr matrix
#' @param mar margin to apply the transformation (rows: 1 or cols: 2)
#' @rdname alr
#' @export
alr.matrix <- function(x.f, mar=2, divcomp=1, base=exp(1), removeDivComp=TRUE,
                        tol=.Machine$double.eps, ...) {
    if (mar == 1) x.f <- t(x.f)
    zero <- (x.f >= tol)
    LOG <- log(ifelse(zero, x.f, 1), base)
    x.alr <- ifelse(zero, LOG - LOG[,divcomp], 0.0)
    if (removeDivComp) x.alr[,-divcomp]
    else x.alr
}

#' @method alr data.frame
#' @rdname alr
#' @export
alr.data.frame <- function(x.f, mar=2, ...) {
    alr(as.matrix(x.f), mar, ...)
}

#' @keywords internal
triu <- function(x) x[upper.tri(x)]
#' @keywords internal
tril <- function(x) x[lower.tri(x)]

#' @keywords internal
triu2diag <- function(x, diagval=0) {
    e <- length(x)
    n <- .5 * (sqrt(8*e + 1)+1)
    mat <- matrix(0, n, n)
    mat[upper.tri(mat)] <- x
    mat <- mat + t(mat)
    diag(mat) <- diagval
    mat
}
