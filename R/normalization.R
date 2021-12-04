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


#' Normalize via dirichlet sampling
#'
#' "Normalize" a count vector by drawing a single sample from a Dirichlet distribution, using the count vector as the prior.
#' @param x count data vector
#' @export
norm_rdiric <- function(x) {
  VGAM::rdiric(n=1, shape=x)[1,]
}

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

#' Centered log-ratio functions
#' @param x input count data
#' @param pseudo value of the pseudocount
#' @param ... pass through arguments
#' @rdname clr
#' @export
pclr <- function(x, mar=2, pseudo=1, ...) {
    clr(x+pseudo, mar=mar, ...)
}


#' Modified central log ratio (mclr) transformation
#'
#' Introduced by GraceYoon/SPRING
#'
#' @param data raw count data or compositional data (n by p) does not matter.
#' @param base exp(1) for natural log
#' @param tol tolerance for checking zeros

# For eps and atleast, users do not have to specify any values. Default should be enough.
#' @param eps epsilon in eq (2) of the paper "Yoon, Gaynanova, M\"{u}ller (2019), Frontiers in Genetics". positive shifts to all non-zero compositions. Refer to the paper for more details. eps = absolute value of minimum of log ratio counts plus c.
#' @param atleast default value is 1. Constant c which ensures all nonzero values to be strictly positive. default is 1.
#'
#'
#' @return \code{mclr} returns a data matrix of the same dimension with input data matrix.
#' @export
## from GraceYoon/SPRING/src/R/helpers.R
mclr <- function(data, base = exp(1), tol = 1e-16, eps = NULL, atleast = 1){
  data <- as.matrix(data)
  nzero <- (data >= tol)  # index for nonzero part
  LOG <- ifelse(nzero, log(data, base), 0.0) # take log for only nonzero values. zeros stay as zeros.

  # centralize by the log of "geometric mean of only nonzero part" # it should be calculated by each row.
  if (nrow(data) > 1){
    clrdata <- ifelse(nzero, LOG - rowMeans(LOG)/rowMeans(nzero), 0.0)
  } else if (nrow(data) == 1){
    clrdata <- ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
  }

  if (is.null(eps)){
    if(atleast < 0){
      warning("atleast should be positive. The functions uses default value 1 instead.")
      atleast = 1
    }
    if( min(clrdata) < 0 ){ # to find the smallest negative value and add 1 to shift all dataa larger than zero.
      positivecst <- abs(min(clrdata)) + atleast # "atleast" has default 1.
    }else{
      positivecst <- 0
    }
    # positive shift
    ADDpos <- ifelse(nzero, clrdata + positivecst, 0.0) ## make all non-zero values strictly positive.
    return(ADDpos)
  } else if(eps == 0){
    ## no shift. clr transform applied to non-zero proportions only. without pseudo count.
    return(clrdata)
  } else if(eps > 0){
    ## use user-defined eps for additional positive shift.
    ADDpos <- ifelse(nzero, clrdata + eps, 0.0)
    return(ADDpos)
  } else {
    stop("check your eps value for additional positive shift. Otherwise, leave it as NULL.")
  }
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

# TODO: should these be exported
#' @keywords internal
dclr <- function(x) t(clr(apply(x, 1, norm_diric),2))

#' @keywords internal
dclrNPN <- function(x) huge::huge.npn(t(clr(apply(x, 1, norm_diric),2)), verbose=FALSE)


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
