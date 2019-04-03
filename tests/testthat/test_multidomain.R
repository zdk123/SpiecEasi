context('test multidomain spiec.easi')


p <- 10
e <- p
n <- 100
set.seed(10010)
g <- make_graph('erdos_renyi', p, e)
S <- cov2cor(prec2cov(graph2prec(g)))
X <- exp(rmvnorm(n, rep(0,p), S))


test_that("execution succeeds", {
  expect_s3_class(
    suppressWarnings(
      spiec.easi(list(X[,1:5], X[,6:10]), method='mb', nlambda=4,
        verbose=FALSE, lambda.min.ratio=1e-1, pulsar.select=FALSE)),
    "pulsar.refit")
})
