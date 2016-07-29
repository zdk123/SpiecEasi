#' get opt merge
#'
#' Get the optimal merge matrix StARS is run.
#'
#' @param est output from \code{spiec.easi} or \code{icov.select}
#' @return Weighted, symmetric matrix of edge probability estimates.
#' @export
getOptMerge <- function(est) {
    if (class(est) == "select" && est$criterion == "stars") {
        return(est$merge[[est$opt.index]])
    } else 
        stop("Run spiec-easi with criterion=\"stars\"")

}


#' get opt beta
#'
#' Get the optimal beta matrix when MB neighborhood selection is run with StARS
#'
#' @param est output from \code{spiec.easi} or \code{icov.select}
#' @return Weighted, non-symmetric matrix of model coefficients.
#' @export
getOptBeta <- function(est) {
    if (class(est) == "select" && est$method == "mb") {
        return(est$beta[[est$opt.index]])
    } else 
        stop("Run spiec-easi with method=\"mb\"")
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
#'  \item ave:    Arithmetic average between the two possible values of beta
#'  \item maxabs: The maximum [absolute] value between the two values
#'  \item upper:  Take the values from the upper triangle
#'  \item lower:  Take the values from the lower triangle
#'}
#' @return a symmetric coefficient matrix
#' @export
symBeta <- function(beta, mode='ave') {
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
