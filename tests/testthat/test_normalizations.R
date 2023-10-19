context('setup')
## Test various correlation methods
p <- 20
e <- p
n <- 100
set.seed(10010)
g <- make_graph('erdos_renyi', p, e)
S <- cov2cor(prec2cov(graph2prec(g)))
X <- rmvnegbin(n, rep(10,p), S, ks=.3)

test_that("execution succeeds", {
  Xmclr <- .spiec.easi.norm(X, method='mclr')
  Xpclr <- .spiec.easi.norm(X, method='pclr')
  Xalr <- .spiec.easi.norm(X, method='alr')
})
