#' Adjacency to igraph
#'
#' Convert an adjacency matrix (ie - from the \code{sparseiCov} function) to an igraph object
#'
#' @param Adj an Adjacency matrix
#' @param rmEmptyNodes should unconnected nodes be removed from the graph
#' @param diag Flag to include self-loops (diagonal of adjacency matrix)
#' @param edge.attr named list of attributes for graph edges
#' @param vertex.attr named list of attributes for graph vertices
#' @return An igraph object
#' @export
#' @examples
#' # Create a symmetric adjacency matrix
#' adj <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow=3, byrow=TRUE)
#' 
#' # Convert to igraph
#' g <- adj2igraph(adj, vertex.attr=list(name=c('A', 'B', 'C')))
adj2igraph <- function(Adj, rmEmptyNodes=FALSE, diag=FALSE, edge.attr=list(),
                       vertex.attr=list(name=1:ncol(Adj))) {
    g <- igraph::graph_from_adjacency_matrix(Adj, mode = "undirected", weighted = TRUE, diag=diag)

    if (length(vertex.attr) > 0) {
        for (i in seq_along(vertex.attr)) {
            attr <- names(vertex.attr)[i]
            g <- igraph::set_vertex_attr(g, attr, index=igraph::V(g), vertex.attr[[i]])
        }
    }

    if (length(edge.attr) > 0) {
        for (i in seq_along(edge.attr)) {
            attr <- names(edge.attr)[i]
            g <- igraph::set_edge_attr(g, attr, index=igraph::E(g), edge.attr[[i]])
        }
    }

    if (rmEmptyNodes) {
        ind <- igraph::V(g)$name[which(igraph::degree(g) < 1)]
        g <- igraph::delete_vertices(g, ind)
#        names <- names[-ind]
    }
    g
}
