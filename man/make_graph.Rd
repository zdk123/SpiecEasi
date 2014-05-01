\name{make_graph}
\alias{make_graph}
\title{Procedure to generate graph topologies for Gaussian Graphical Models}
\usage{
make_graph(method, D, e, enforce = TRUE, ...)
}
\arguments{
  \item{method}{Type of graph to make}

  \item{D}{Number of nodes/OTUs (Graph dimension)}

  \item{e}{Number of edges (preferably sparse, must be at
  least 1/2 D)}

  \item{enforce}{add/remove edges to enforce graph has e
  edges}

  \item{...}{additional options to graph method}
}
\description{
Procedure to generate graph topologies for Gaussian
Graphical Models
}

