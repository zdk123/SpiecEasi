context('setup')

p <- 20
e <- p
n <- 500
set.seed(10010)
g <- make_graph('erdos_renyi', p, e)
S <- cov2cor(prec2cov(graph2prec(g)))
X <- exp(rmvnorm(n, rep(0,p), S))

pargs <- list(seed=10010, rep.num=10)

context("SPIEC-EASI fit")
lmx  <- .7
lmr  <- 1e-2
nlam <- 20
## No selection ##
out <- spiec.easi(X, method='mb', verbose=FALSE, lambda.max=lmx, lambda.min.ratio=lmr, nlambda=nlam, pulsar.select=FALSE)

## StARs / B-StARS
t1 <- system.time(
out.stars <- spiec.easi(X, method='mb', verbose=FALSE, lambda.max=lmx,
  lambda.min.ratio=lmr, nlambda=nlam, sel.criterion='stars',
  pulsar.select=TRUE, pulsar.params=pargs))
t2 <- system.time(
out.bstars <- spiec.easi(X, method='mb', verbose=FALSE, lambda.max=lmx,
  lambda.min.ratio=lmr, nlambda=nlam, sel.criterion='bstars',
  pulsar.select=TRUE, pulsar.params=pargs))
## Batch Mode StARs / B-StARS
options(batchtools.verbose=FALSE)
bout.stars <- spiec.easi(X, method='mb', verbose=FALSE, lambda.max=lmx,
  lambda.min.ratio=lmr, nlambda=nlam, sel.criterion='stars',
  pulsar.select='batch', pulsar.params=pargs)
bout.bstars <- spiec.easi(X, method='mb', verbose=FALSE, lambda.max=lmx,
  lambda.min.ratio=lmr, nlambda=nlam, sel.criterion='bstars',
  pulsar.select='batch', pulsar.params=pargs)


test_that("no pulsar has same output", {
  tmp <- out.stars$select$stars$opt.index
  expect_equal(as.matrix(Matrix::drop0(out$est$path[[tmp]])),
               as.matrix(out.stars$refit$stars))
})

test_that("stars == bstars", {
  expect_equal(as.matrix(out.bstars$refit$stars),
               as.matrix(out.stars$refit$stars))
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

context('spiec.easi getters')

test_that("Getter API throws errors if no pulsar selection", {
  expect_error(getOptInd(out))
  expect_error(getOptNet(out))
  expect_error(getRefit(out))
  expect_error(getOptLambda(out))
  expect_error(getOptMerge(out))
  expect_error(getOptCov(out))
  expect_error(getOptiCov(out))
  expect_error(getOptBeta(out))
})

runtests <- function(out) {
  expect_equal(getOptInd(out),
              (i<-out$select$stars$opt.index))
  expect_equal(getOptNet(out),
               Matrix::drop0(out$refit$stars))
  expect_equal(getRefit(out),
               Matrix::drop0(out$refit$stars))
  expect_equal(getOptLambda(out), out$lambda[i])
  expect_equal(getOptMerge(out),
               Matrix::drop0(out$select$stars$merge[[i]]))
  expect_equal(getOptBeta(out),
               Matrix::drop0(out$est$beta[[i]]))
}

test_that("Getter API, pulsar / stars ", {
  runtests(out.stars)
})

test_that("Getter API, pulsar / bstars ", {
  runtests(out.bstars)
})

test_that("Getter API, batch pulsar / stars ", {
  runtests(bout.stars)
})

test_that("Getter API, batch pulsar / bstars ", {
  runtests(bout.bstars)
})
