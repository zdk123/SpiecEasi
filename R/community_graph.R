#' Precision matrix (inverse covariance) to a covariance matrix
#'
#' @param Precision symmetric precision matrix
#' @param tol tolerance to define a zero eigenvalue (ie - is Prec positive definite)
#' @importFrom MASS ginv
#' @export
prec2cov <- function(Precision, tol=1e-4) {
    eigval <- eigen(Precision)$values
    if (any(eigval < tol)) {
        warning("Warning: Precision matrix not invertible, trying generalized inverse instead")
        ginv(Precision)
    } else {
        solve(Precision)
    }
}

#' Covariance matrix to its matrix inverse (Precision matrix)
#'
#' @param Cov symmetric covariance matrix (can be correlation also)
#' @param tol tolerance to define a zero eigenvalue (ie - is Prec positive definite)
#' @importFrom MASS ginv
#' @export
cov2prec <- function(Cov, tol=1e-4) {
    Precision <- tryCatch(solve(Cov), error = function(e) {
              warning("Warning: Precision matrix not invertible, trying generalized inverse instead")
              ginv(Cov) })
    Precision[which(Precision == 0)] <- tol
}




#' Convert a symmetric graph (extension of R matrix class)
#'
#'  Has internal rules for converting various graph topologies into the associated
#' adjancency and, therefore, precision matrix
#'
#' @export
#' @param Graph graph adjacency matrix
#' @param posThetaLims length 2 vector of lower and upper bound of positive values
#' @param negThetaLims length 2 vector of lower and upper bound of negative values
#' @param targetCondition sets the condition of the precision matrix by modulating the magnitude of the diagonal
#' @param epsBin the convergence tolerance of the condition number binary search
#' @param numBinSearch maximum number of iterations
graph2prec <- function(Graph, posThetaLims=c(2,3), negThetaLims=-posThetaLims, targetCondition=100, epsBin=1e-2,
                        numBinSearch=100) {

    if (class(Graph) != 'graph') stop('input is not a graph')
    n <- ncol(Graph)

    posThetaLims <- sort(posThetaLims)
    negThetaLims <- sort(negThetaLims)

    #add diagonal
    Theta <- Graph + diag(n)
    # node degree
    degVec <- colSums(Graph)

    # add random edge weights uniformly from theta_min to theta_max
    utri  <- triu(Theta)
    nzind <- which(utri != 0)

    if (length(posThetaLims) > 2 || length(negThetaLims) > 2)
        stop("theta_max and theta_min should be a numeric vector of length 1 or 2")

    rands <- runif(length(nzind))
    mapToRange <- function(x, lim) {
        span <- diff(sort(lim))
        min(lim) + (x * span)
    }

    boolind <- sample(c(TRUE, FALSE), length(nzind), replace=TRUE)

    rands[boolind]  <- mapToRange(rands[boolind], posThetaLims)
    rands[!boolind] <- mapToRange(rands[!boolind], negThetaLims)

    utri[nzind] <- rands
    Theta <- triu2diag(utri, 1)

    # find smallest eigenvalue such that Theta is invertible
    eigVals <- eigen(Theta)$values
    minEig  <- min(eigVals)
    maxEig  <- max(eigVals)

    if (minEig < 1e-2) Theta <- Theta + abs(minEig)*diag(n)
    diagConst <- .binSearchCond(Theta, targetCondition, numBinSearch, epsBin)
    Theta <- Theta + diagConst*diag(n)
    return(Theta)
}


.binSearchCond <- function(Theta, condTheta, numBinSearch, epsBin) {
# Internal function that determines the constant in the diagonal to satisfy the
# condition constraint on the Precision/Covariance matrix

    n <- nrow(Theta)
    currCondTheta <- kappa(Theta)
    if (currCondTheta < condTheta) {
        # Max entry in the diagonal (lower bound)
        currLB   <- -max(diag(Theta))
        stepSize <- currLB+.Machine$double.eps

        while (currCondTheta < condTheta) {
            currCondTheta <- kappa(Theta+stepSize*diag(n))
            stepSize      <- stepSize/2
        }
        currUB <- stepSize
    } else {
        currLB <- 0
        stepSize = 0.1

        while (currCondTheta > condTheta) {
            currCondTheta <- kappa(Theta + stepSize*diag(n))
            stepSize      <- 2*stepSize
        }
        currUB <- stepSize
    }

    for (i in 1:numBinSearch) {
        diagConst <- (currUB+currLB)/2
        currCondTheta <- kappa(Theta+diagConst*diag(n))

        if (currCondTheta < condTheta) currUB <- diagConst
        else currLB <- diagConst

        if (abs(currCondTheta-condTheta)<epsBin) break
    }
    diagConst
}


#' Procedure to generate graph topologies for Gaussian Graphical Models
#' @param method Type of graph to make
#' @param D Number of nodes/OTUs (Graph dimension)
#' @param e Number of edges (preferably sparse, must be at least 1/2 D)
#' @param enforce add/remove edges to enforce graph has e edges
#' @param ... additional options to graph method
#' @export
make_graph <- function(method, D, e, enforce=TRUE, ...) {
    if (e>((D-1)*D)/2) stop('Number of edges e must smaller than D(D-1)/2')

    method <- switch(method, cluster = "cluster", erdos_renyi = "erdos_renyi",
                       hub = "hub", scale_free = "scale_free",
                       block = "block", band = "band",
                       stop(paste("Error: graph method", method, "not supported")))

    if (method != "erdos_renyi")
      if (e < round(D/2)) stop('Number of edges e must be bigger than 1/2 D')
    graphgen <- get(method)
    Graph    <- graphgen(D, e=e, ...)
    attr(Graph, "graph") <- method

    ## enforce edge number
    if (enforce) {
        Graph <- enforceE(Graph, e)
    }
    return(structure(Graph, class='graph'))
}


#' @keywords internal
enforceE <- function(Graph, e) {

    nedges <- edge_count(Graph)
    D      <- nrow(Graph)
    diffE  <- nedges - e
    uniqMatInds <- function(inds, D) {
        tmpInds <- arrayInd(inds, c(D,D))
        unique(t(unique(apply(tmpInds,1,sort))))
    }
    if (diffE > 0) {
        oneInds <- which(Graph == 1)
        tmpInds <- uniqMatInds(oneInds, D)
        randi   <- sample(1:nrow(tmpInds), abs(diffE))
        for (i in randi) {
            gind <- tmpInds[i,]
            Graph[gind[1], gind[2]] <- Graph[gind[2], gind[1]] <- 0

        }
    } else if (diffE < 0) {
        diag(Graph) <- 1   # fill diag with dummy 1s
        zeroInds <- which(Graph == 0)
        tmpInds <- uniqMatInds(zeroInds, D)
        randi   <- sample(1:nrow(tmpInds), abs(diffE))
        for (i in randi) {
            gind <- tmpInds[i,]
            Graph[gind[1], gind[2]] <- Graph[gind[2], gind[1]] <- 1
        }
        diag(Graph) <- 0
    }

    if (diffE != 0) enforceE(Graph, e)
    else Graph
}


#' @keywords internal
scale_free <- function(D, e, pfun) {
# Make a scale free graph
# Args:
#   D -   > The number of nodes
#   pfun -> override default probability of getting an edge

    if (e < (D-1)) stop("Too few edges to generate a scale-free graph, e>=D")

    #initialize D by D zero matrix
    Graph <- matrix(0, D, D)
    K_mat <- matrix(0, D)   # keep track of the number of degrees in the graph

    # pick first two nodes at random
    nodes <- sample(1:D, 2)
    #connect nodes
    Graph[nodes[1], nodes[2]] <- 1
    Graph[nodes[2], nodes[1]] <- 1
    # Add degree to nodes
    K_mat[nodes[1]] <- K_mat[nodes[1]] + 1
    K_mat[nodes[2]] <- K_mat[nodes[2]] + 1

    if (missing(pfun)) {
        pfun <- function(K_mat, i) {
            (K_mat[i] / sum(K_mat))
        }
    }

    newnodes <- setdiff(1:D, nodes)
    existnodes <- nodes
    while (length(newnodes) != 0) {  # For each node
        n1 <- na.exclude(newnodes)[1]
        for (n2 in existnodes)  { # and a node already in the graph
            p <- pfun(K_mat, n2)
            conn <- sample(c(TRUE, FALSE), 1, prob=c(p, 1-p))
            if (conn) {
                 #connect nodes
                Graph[n1, n2] <- 1
                Graph[n2, n1] <- 1
                 # Add degree to nodes
                K_mat[n1] <- K_mat[n1] + 1
                K_mat[n2] <- K_mat[n2] + 1
                # add connected nodes to list, remove from new
                existnodes <- c(existnodes, n1)
                newnodes   <- setdiff(newnodes, n1)
                break
            } else {
                # move n1 to back of the list
                if (length(newnodes > 1)) newnodes <- newnodes[c(2:length(newnodes),1)]
            }
        }
    }
    return(Graph)
}

#' @keywords internal
erdos_renyi <- function(D, e, p=e/(D*(D-1)/2)) {
# Make a random graph:
# Args:
#   D -> number of nodes
#   e -> Number of edges
#   p -> probability of edge
    if (missing(e)) stop("Must supply either number of edges (e)")
    Gtemp <- matrix(runif(D^2, 0,1), D, D)
    Gtemp <- ifelse(Gtemp < p, 1, 0)
    Gtemp[lower.tri(Gtemp, diag=TRUE)] <- 0
    Graph <- Gtemp + t(Gtemp)
    return(Graph)
}

#' @keywords internal
hub <- function(D, e, numHubs=ceiling(D/20)) {
# Make hub graph
# Args:
# D   -> number of nodes
# e   -> number of edges, only for compatability
# numHubs -> number of hubs (divide nodes into this many groups)

    Graph   <- matrix(0, D, D)
    # partiion nodes into groups
    groupid  <- factor(sample(1:numHubs, D, replace=TRUE))
    groupind <- split(1:D, groupid)

    for (ind in groupind) {
        hub <- sample(ind, 1)   # pick 1 hub
        nodes <- setdiff(ind, hub)
        Graph[hub, nodes] <- 1
        Graph[nodes, hub] <- 1
    }
    attr(Graph, "g") <- numHubs
    return(Graph)
}

#' @keywords internal
cluster <- function(D, e, numHubs=floor((D/15)+(e/D))-1) {
# Make a cluster graph (groups of random graphs)
# Args:
# D   -> number of nodes
# numHubs -> number of hubs (divide nodes into this many groups)
# d    -> number of edges
    Graph   <- matrix(0, D, D)
    # partiion nodes into groups
#    groupid  <- factor(sample(1:numHubs, D, replace=TRUE))
    groupid <- factor(sample(rep(1:numHubs, length.out=D)))
    groupind <- split(1:D, groupid)
    eHub     <- ceiling(e/numHubs)
    for (ind in groupind) {
        nHub <- length(ind)
        p    <- eHub/(nHub*(nHub-1)/2)
        Graph[ind,ind] <- erdos_renyi(length(ind), e=eHub, p=p)
    }
    attr(Graph, "g") <- numHubs
    return(Graph)
}

#' @keywords internal
band <- function(D, e) {
# Make a banded graph
# Args:
#   D    -> number of nodes
#   e    -> target number of edges

    bestFit <- FALSE
    k <- 1
    Graph   <- matrix(0, D, D)
    zero    <- matrix(0, D, D)
    # add off diagonals until e is exhausted
    while (!bestFit) {
        off1    <- matrix(0, D, D)
        diag(off1[-(1:k),]) <- rep(1, D-k)
        off2    <- matrix(0, D, D)
        diag(off2[,-(1:k)]) <- rep(1, D-k)
        tempG <- Graph + off1 + off2
        if (sum(tempG)>(2*e)) bestFit = TRUE
        else Graph <- tempG
        k <- k + 1
    }
    return(Graph)
}

#' @keywords internal
block <- function(D, e, numHubs) {  #blocksize=20, p=D/((D/blocksize)*(blocksize*(blocksize)/2)), u=NULL, v=NULL) {
# Make precision matrix for block graph (note the difference from other functions)
# Args:
#   D    -> Number of nodes
#   e    -> target number of edges
#   numHubs -> number of hubs, calculate best number if missing

    if (missing(numHubs)) {
        bestFit <- Inf
        for (numHubs in 2:(D/2)) {
            nHubs   <- round(D/numHubs)
            currFit <- abs(nHubs*(nHubs-1)/2*numHubs - e)
            if (currFit < bestFit) {
                bestNumHubs <- numHubs
                bestNHubs   <- nHubs
                bestFit     <- currFit
            }
        }
        numHubs <- bestNumHubs
    }

    Graph <- matrix(0, D, D)
   # partition nodes into blocks
    blockid  <- factor(sample(1:numHubs, D, replace=TRUE))
    blockind <- split(1:D, blockid)
  #  hubs <- round(seq.int(1,D+1, length.out=numHubs+1))
    for (ind in blockind) {
        nHubs <- length(ind)
        Ghub  <- matrix(1, nHubs, nHubs)
        Ghub[seq(1, nHubs^2, nHubs+1)] <- 0
        Graph[ind,ind] <- Ghub
    }
    return(Graph)
}


#' @keywords internal
edge_count <- function(Graph) {
    # return number of edges in a [square matrix] graph
    length(which(Graph[upper.tri(Graph)] != 0))
}

#' @keywords internal
eigGraph <- function(G, tol=1e-12) {
    D <- diag(rowSums(G))
    L <- D-G
    eigsG   <- eigen(L)$values
    specGap <- eigsG[which(eigsG > tol)]
    nCC <- sum(eigsG < tol)
    return(list(specGap=specGap, nCC=nCC, eigsG=eigsG))
}

#' @keywords internal
graphReport <- function(Graph) {
    nedges   <- edge_count(Graph)
    eigGraph <- eigGraph(Graph)
    return(c(list(NumEdges=nedges), eigGraph))
}

#' @keywords internal
covReport <- function(Cov, Prec) {
    condCov <- kappa(Cov)
    condPrec <- kappa(Prec)
    # Stats over correlations
    corr  <- cov2cor(Cov)
    corr <- corr[corr <0.999]
    corr <- corr[corr != 0]
    minC  <- min(corr)
    meanC <- mean(corr)
    medianC <- median(corr)
    maxC  <- max(corr)
    corrStats <- t(data.frame(minC, meanC, medianC, maxC))
    return(list(condCov=condCov, condInvCov=condPrec, corrStats=corrStats))
}


#' s3 method for graph to other data types
#' @param x graph adjacency matrix
#' @param ... Arguments to base as.data.frame
#' @export
as.data.frame.graph <- function(x, ...) {
    as.data.frame(as.matrix(x), ...)
}

#' s3 method for graph to other data types
#' @param x graph adjacency matrix
#' @param ... Arguments to base as.matrix
#' @export
as.matrix.graph <- function(x, ...) {
    class(x) <- 'matrix'
    return(x)
}
