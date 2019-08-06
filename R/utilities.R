#' get StARS-optimal network
#'
#' Get the optimal network, and related structures, when StARS is run.
#'
#' @param est output from \code{spiec.easi}
#' @return numeric or matrix associated with a StARS solution.
#' @details
#' Use the getter functions to parse \code{\link{spiec.easi}} output:
#'
#' \itemize{
#'   \item getOptLambda: penalty parameter from provided lambda path
#'   \item getOptInd: index of the selected lambda from provided lambda path
#'   \item getOptNet / getRefit: the optimal (StARS-refit) network
#'   \item getStability: average stability at the selected sparsity
#'   \item getOptMerge: symmetric matrix with edge-wise stability
#'   \item getOptiCov: the optimal inverse covariance matrix (glasso only)
#'   \item getOptCov: the optimal covariance matrix associated with the selected network (glasso only)
#'   \item getOptBeta: the optimal coefficient matrix (mb only)
#'}
#' @export
getOptInd <- function(est)   getOptX(est, 'index')

#' @rdname getOptInd
#' @export
getOptLambda <- function(est)   getOptX(est, 'lambda')

#' @rdname getOptInd
#' @export
getOptMerge <- function(est) getOptX(est, 'merge')

#' @rdname getOptInd
#' @export
getStability <- function(est)  getOptX(est, 'stars')

#' @rdname getOptInd
#' @export
getOptNet <- function(est)   getOptX(est, 'refit')

#' @rdname getOptInd
#' @export
getRefit <- function(est)    getOptNet(est)

#' @rdname getOptInd
#' @export
getOptBeta <- function(est)  getOptX(est, 'beta')

#' @rdname getOptInd
#' @export
getOptCov <- function(est)   getOptX(est, 'cov')

#' @rdname getOptInd
#' @export
getOptiCov <- function(est)  getOptX(est, 'icov')

#' @keywords internal
getOptX <- function(est, ...) {
  UseMethod('getOptX', est)
}

#' @keywords internal
getOptX.pulsar.refit <- function(est, getter='index') {
  if (length(est$refit) > 0) {
    if (names(est$refit) %in% "stars") {
      return(switch(getter,
        index = {
              est$select$stars$opt.index
                },
        refit = {
              Matrix::drop0(est$refit$stars)
        },
        merge = {
              Matrix::drop0(est$select$stars$merge[[getOptInd(est)]])
        },
        stars = {
              est$select$stars$summary[getOptInd(est)]
        },
        lambda = {
              est$lambda[[getOptInd(est)]]
        },
        icov = {
          if (est$est$method == "glasso")
            Matrix::drop0(est$est$icov[[getOptInd(est)]])
          else
          stop("Run spiec-easi with method=\"glasso\"")
        },
        cov = {
          if (est$est$method == "glasso")
            Matrix::drop0(est$est$cov[[getOptInd(est)]])
          else
            stop("Run spiec-easi with method=\"glasso\"")
        },
        beta = {
          if (est$est$method == "mb")
            Matrix::drop0(est$est$beta[[getOptInd(est)]])
          else
            stop("Run spiec-easi with method=\"mb\"")
        }
      ))
    } else
      stop("Run spiec.easi with sel.criterion=\"stars\"")
  } else
    stop("Run spiec.easi with pulsar.select=TRUE")
}



#' sym beta
#'
#' Symmetrize a beta (coefficient) matrix, ie. selected from MB neighborhood selection
#'
#' @param beta square coefficient matrix
#' @param mode how to symmetrize, see details
#' @details
#' Mode can be:
#' \itemize{
#'  \item ave:  Arithmetic average between the two possible values of beta
#'  \item maxabs: The maximum [absolute] value between the two values
#'  \item upper:  Take the values from the upper triangle
#'  \item lower:  Take the values from the lower triangle
#'}
#' @return a symmetric coefficient matrix
#' @export
symBeta <- function(beta, mode='ave') {
  t <- Matrix::t
  if (nrow(beta) != ncol(beta)) {
    stop('expecting a symmetric Matrix')
  }
  if (mode=='ave') {
    symbeta <- (beta+t(beta))/2
   } else if (mode == "maxabs") {
    upt <- Matrix::triu(beta)
    lot <- t(Matrix::tril(beta))
    suppressMessages(maxt <- pmax(abs(upt), abs(lot)))
    uptind <- Matrix::which(maxt == abs(upt))
    lotind <- Matrix::which(maxt == abs(lot))
    if (length(uptind != 0)) maxt[uptind] <- maxt[uptind]*sign(upt[uptind])
    if (length(lotind != 0)) maxt[lotind] <- maxt[lotind]*sign(lot[lotind])
    symbeta <-  maxt + t(maxt)
  } else if (mode == "upper") {
    upt <- Matrix::triu(beta)
    symbeta <-  upt + t(upt)
  } else if (mode == "lower") {
    lot <- Matrix::tril(beta)
    symbeta <- lot + t(lot)
  } else
    stop ("mode not recognized")
  as(symbeta, 'symmetricMatrix')
}
