#' stars.roc, stars.pr
#'
#' Plot a ROC (reciever operator characteristic) or a Precision-Recall curve
#' along the stars 'confidence path'. Each edge is a number in [0,1],
#' which is on the fraction of inferred graphs over subsamples in which that edge
#' appeared in stars.
#'
#' @param optmerge the optimal 'merge' matrix selected by stars
#' @param theta the true graph or precision matrix
#' @param verbose display messages
#' @param plot graph the output
#' @param ll number of points for the plot
#' @importFrom grDevices dev.off png
#' @export
stars.roc <- function(optmerge, theta, verbose = TRUE, plot = TRUE, ll=15) {
    if (!plot) { ff <- tempfile() ; png(filename=ff) }
    tmp <- huge::huge.roc(merge2path(optmerge, ll), theta, verbose)
    if (!plot) { dev.off() ; unlink(ff) }
    return(tmp)
}

#' @rdname stars.roc
#' @export
stars.pr <- function(optmerge, theta, verbose = TRUE, plot = TRUE, ll=15) {
    huge.pr(merge2path(optmerge, ll+1), theta, verbose, plot)
}

merge2path <- function(merge, length.out) {
 if (missing(length.out)) {
    path <- unique(c(merge))
 } else {
    path <- seq(min(merge), max(merge), length.out=length.out)
 }
 lapply(path, function(i) merge > i)
}

#' @noRd
huge.pr <- function (path, theta, verbose = TRUE, plot = TRUE) {
    gcinfo(verbose = FALSE)
    ROC = list()
    theta = as.matrix(theta)
    d = ncol(theta)
    pos.total = sum(theta != 0)
    neg.total = d * (d - 1) - pos.total
    if (verbose)
        message("Computing F1 scores, false positive rates and true positive rates....", appendLF=FALSE)
    ROC$prec = rep(0, length(path))
    ROC$rec  = rep(0, length(path))
    ROC$F1 = rep(0, length(path))
    for (r in 1:length(path)) {
        tmp = as.matrix(path[[r]])
        tp.all = (theta != 0) * (tmp != 0)
        diag(tp.all) = 0
        ROC$tp[r] <- sum(tp.all != 0)/pos.total
        fp.all = (theta == 0) * (tmp != 0)
        diag(fp.all) = 0
        ROC$fp[r] <- sum(fp.all != 0)/neg.total
        fn = 1 - ROC$tp[r]
        precision = ROC$tp[r]/(ROC$tp[r] + ROC$fp[r])
        recall = ROC$tp[r]/(ROC$tp[r] + fn)
        ROC$prec[r] <- precision
        ROC$rec[r]  <- recall
        ROC$F1[r] = 2 * precision * recall/(precision + recall)
        if (is.na(ROC$F1[r]))
            ROC$F1[r] = 0
    }
    if (verbose)
        message("done.")
    rm(precision, recall, tp.all, fp.all, path, theta, fn)
    gc()
    ord.p = order(ROC$prec, ROC$rec, na.last=NA)
    ROC$prec <- ROC$prec[ord.p]
    ROC$rec  <- ROC$rec[ord.p]
    tmp2 = c(min(c(4/(d-1), ROC$prec)), ROC$prec[ord.p], 1)
    tmp1 = c(1, ROC$rec[ord.p], 0)
    if (plot) {
        par(mfrow = c(1, 1))
        plot(tmp1, tmp2, type = "b", main = "PR Curve", xlab = "Recall",
            ylab = "Precision", ylim = c(0, 1))
    }
    tmax <- diff(range(tmp2))*diff(range(tmp1))
    ROC$AUC = sum(diff(tmp2) * (tmp1[-1] + tmp1[-length(tmp1)]))/(2*tmax)
    rm(ord.p, tmp1, tmp2)
    gc()
    class(ROC) = "roc"
    return(ROC)
}

#' Edge set dissimilarity
#'
#' Compute the dissimilarity between the edge sets of two networks via:
#' \enumerate{
#'  \item maximum overlap: |x  y| / max\{|x|,|y|\}
#'  \item jaccard index (default):   |x  y|/(|x U y|)
#' }
#' Input networks do not have to have the same node sets.
#' @param x pxp adjacency matrix (\code{Matrix::sparseMatrix} class)
#' @param y other qxq adjacency matrix (\code{Matrix::sparseMatrix} class)
#' @param metric 'jaccard' or 'max'
#' @param otux taxa names of adjacency x
#' @param otuy taxa names of adjacency y
#' @export
edge.diss <- function(x, y, metric='jaccard', otux=NULL, otuy=NULL) {
  stopifnot(inherits(x, 'sparseMatrix'), TRUE)
  stopifnot(inherits(y, 'sparseMatrix'), TRUE)
  xli <- Matrix::summary(Matrix::triu(x))
  yli <- Matrix::summary(Matrix::triu(y))
  if (!is.null(otux)) {
    xli[,1] <- otux[xli[,1]]
    xli[,2] <- otux[xli[,2]]
  }
  if (!is.null(otuy)) {
    yli[,1] <- otuy[yli[,1]]
    yli[,2] <- otuy[yli[,2]]
  }
  xedges <- apply(xli[,1:2], 1, paste, collapse="-")
  yedges <- apply(yli[,1:2], 1, paste, collapse="-")
  if (metric=="jaccard") {
    return(length(intersect(xedges, yedges)) / length(unique(c(xedges, yedges))))
  } else {
    return(length(intersect(xedges, yedges)) / max(length(xedges), length(yedges)))
  }
}