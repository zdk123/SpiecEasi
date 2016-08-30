library(profvis)

getNet <- function(type, p, e) {
    G <- diag(p) - 1/p
    method <- type
    graph <- enforceE(match.fun(method)(p, e), e)
    class(graph) <- "graph"
    if (method == "Null") { Prec <- as.matrix(graph) ; diag(Prec) <- 3 }
    else Prec  <- graph2prec(graph)
    Cov   <- prec2cov(Prec)
    S     <- G%*%Cov%*%G
    return(list(graph=graph, Prec=Prec, Cov=Cov, S=S))
}

d <- 100
dens <- .05
e <- round(dens*((d*(d-1))/2))
out <- getNet('cluster', d, e)
S <- out$S
graph <- out$graph
Prec <- out$Prec
Cov <- out$Cov

p <- profvis({
  est  <- sparseLowRankiCov(S, r=2, tol=1e-3)
})

print(p)
