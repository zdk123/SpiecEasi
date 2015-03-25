stars.roc <- function(optmerge, theta, verbose = TRUE, ll) {
    huge::huge.roc(merge2path(optmerge, ll), theta, verbose)
}

stars.pr <- function(optmerge, theta, verbose = TRUE, ll) {
    huge.pr(merge2path(optmerge, ll+1), theta, verbose)
}

merge2path <- function(merge, length.out) {

 if (missing(length.out)) {
    path <- unique(c(merge))
 } else {
    path <- seq(min(merge), max(merge), length.out=length.out)
 }
 
 lapply(path, function(i) merge > i)
}

huge.pr <- function (path, theta, verbose = TRUE) {
    gcinfo(verbose = FALSE)
    ROC = list()
    theta = as.matrix(theta)
    d = ncol(theta)
    pos.total = sum(theta != 0)
    neg.total = d * (d - 1) - pos.total
    if (verbose) 
        cat("Computing F1 scores, false positive rates and true positive rates....")
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
        cat("done.\n")
    rm(precision, recall, tp.all, fp.all, path, theta, fn)
    gc()
    ord.p = order(ROC$prec, na.last=NA)
    tmp1 = ROC$prec[ord.p]
    tmp2 = ROC$rec[ord.p]
    par(mfrow = c(1, 1))
    plot(tmp1, tmp2, type = "b", main = "ROC Curve", xlab = "Precision", 
        ylab = "Recall", ylim = c(0, 1))
    ROC$AUC = sum(diff(tmp1) * (tmp2[-1] + tmp2[-length(tmp2)]))/2
    rm(ord.p, tmp1, tmp2)
    gc()
    class(ROC) = "roc"
    return(ROC)
}

