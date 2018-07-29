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
    ord.p = order(ROC$prec, na.last=NA)
    tmp1 = ROC$prec[ord.p]
    tmp2 = ROC$rec[ord.p]
    if (plot) {
        par(mfrow = c(1, 1))
        plot(tmp1, tmp2, type = "b", main = "PR Curve", xlab = "Precision",
            ylab = "Recall", ylim = c(0, 1))
    }
    ROC$AUC = sum(diff(tmp1) * (tmp2[-1] + tmp2[-length(tmp2)]))/2
    rm(ord.p, tmp1, tmp2)
    gc()
    class(ROC) = "roc"
    return(ROC)
}
