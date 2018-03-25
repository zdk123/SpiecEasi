context('setup')

p <- 20
e <- p
n <- 100
set.seed(10010)
g <- make_graph('hub', p, e)
S <- cov2cor(prec2cov(graph2prec(g)))
X <- exp(rmvnorm(n, rep(0,p), S))


pargs <- list(seed=10010, rep.num=10)

context("SPIEC-EASI fit")
lmx  <- .6
lmr  <- 5e-3
nlam <- 35
## No selection ##
out <- spiec.easi(X, method='mb', verbose=FALSE, lambda.max=lmx, lambda.min.ratio=lmr, nlambda=nlam, pulsar.select=FALSE)

## StARs / B-StARS
t1 <- system.time(
out.stars <- spiec.easi(X, method='mb', verbose=FALSE, lambda.max=lmx, lambda.min.ratio=lmr, nlambda=nlam, sel.criterion='stars', pulsar.select=TRUE, pulsar.params=pargs))
t2 <- system.time(
out.bstars <- spiec.easi(X, method='mb', verbose=FALSE, lambda.max=lmx, lambda.min.ratio=lmr, nlambda=nlam, sel.criterion='bstars', pulsar.select=TRUE, pulsar.params=pargs))
## Batch Mode StARs / B-StARS
options(batchtools.verbose=FALSE)
bout.stars <- spiec.easi(X, method='mb', verbose=FALSE, lambda.max=lmx, lambda.min.ratio=lmr, nlambda=nlam, sel.criterion='stars', pulsar.select='batch', pulsar.params=pargs)
bout.bstars <- spiec.easi(X, method='mb', verbose=FALSE, lambda.max=lmx, lambda.min.ratio=lmr, nlambda=nlam, sel.criterion='bstars', pulsar.select='batch', pulsar.params=pargs)


save(out, out.stars, out.bstars, bout.stars, bout.bstars,
    file=system.file('testdata', 'sepulsar.RData', package='SpiecEasi'))

test_that("no pulsar has same output", {
  expect_equal(as.matrix(out$est$path[[out.stars$select$stars$opt.index]]),
               as.matrix(out.stars$refit$stars))
})

test_that("stars == bstars", {
  expect_equal(as.matrix(out.bstars$refit$stars),
               as.matrix(out.stars$refit$stars))
  expect_gt(t1[3], t2[3])

})

test_that("stars == batch stars", {
  expect_equal(as.matrix(bout.stars$refit$stars),
               as.matrix(out.stars$refit$stars))
})

test_that("batch stars == batch bstars", {
  expect_equal(as.matrix(bout.stars$refit$stars),
               as.matrix(bout.bstars$refit$stars))
})

test_that("stars == batch bstars", {
  expect_equal(as.matrix(out.stars$refit$stars),
               as.matrix(bout.bstars$refit$stars))
})


## TODO: test pulsar params in/out
