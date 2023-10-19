context('setup')
## Test various correlation methods
p <- 20
e <- p
n <- 100
set.seed(10010)
g <- make_graph('erdos_renyi', p, e)
S <- cov2cor(prec2cov(graph2prec(g)))
X <- rmvnegbin(n, rep(10,p), S, ks=.3)

pargs <- list(seed=10010, rep.num=10)

context("Correlation fit")
# lmx  <- 1
lmr  <- 5e-2
nlam <- 10

test_that("execution succeeds", {
## Default ##
out_default <- spiec.easi(X, method='mb', cov.method='default', verbose=FALSE, lambda.min.ratio=lmr, nlambda=nlam, pulsar.select=FALSE)
out_cor <- spiec.easi(X, method='mb', cov.method='cor', verbose=FALSE, lambda.min.ratio=lmr, nlambda=nlam, pulsar.select=FALSE)
out_cov <- spiec.easi(X, method='mb', cov.method='cor', verbose=FALSE, lambda.min.ratio=lmr, nlambda=nlam, pulsar.select=FALSE)
out_lcor <- spiec.easi(X, method='mb', cov.method='latentcor', verbose=FALSE, lambda.min.ratio=lmr, nlambda=nlam, pulsar.select=FALSE)

## TODO: test some outputs here
})
