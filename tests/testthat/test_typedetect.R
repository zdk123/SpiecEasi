context('setup')
library(latentcor)


set.seed(10010)
X1 <- gen_data(n = 100, types = c("tru", "con", "con", "bin", "ter"),
                          XP = list(.1,   1e-4, NA,    .5, c(.6,.3,.1)))$X
## create count types
X1[X1[,1]!=0,1] <- round(exp(X1[X1[,1]!=0,1]))
X1[ ,2] <- round(exp(1.5*X1[,2]+3))


p <- 20
e <- p
n <- 100
set.seed(10010)
g <- make_graph('erdos_renyi', p, e)
S <- cov2cor(prec2cov(graph2prec(g)))
X2 <- rmvnegbin(n, rep(10,p), S, ks=.3)

datanorm <- .spiec.easi.norm(list(X1, X2))
types <- attr(datanorm, 'types')

test_that("execution succeeds", {
  est <- spiec.easi(list(X1, X2), method='mb', cov.fun='latentcor')
})
