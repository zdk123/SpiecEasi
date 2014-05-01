#' Convert an adjacency matrix (ie - from the \code{sparseiCov} function) to an igraph object suitable for plotting via phyloseq's \code{plot_network}.
#' 
#' @param Adj an Adjacency matrix
adj2igraph <- function(Adj, names=1:ncol(Adj), rmEmptyNodes=FALSE) {
    g <- graph.adjacency(Adj)

    if (rmEmptyNodes) {
        ind <- which(degree(g) < 1)
        g <- delete.vertices(g, ind)
        names <- names[-ind]
    }

    V(g)$name <- names
    g
}
