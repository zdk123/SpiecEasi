context('setup')

### synth data from American gut

data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- 20 # ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

set.seed(10010)
graph <- make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))


context("marginal attributes preserved")
test_that("Synthetic data has correct attributes: zinegbin", {
    X <- synth_comm_from_counts(amgut1.filt.cs[,1:d], mar=2, distr='zinegbin', Sigma=Cor, n=n)
    expect_equal(ncol(X), d)
    expect_equal(nrow(X), n)
})

test_that("Synthetic data has correct attributes: zipois", {
    X <- synth_comm_from_counts(amgut1.filt.cs[,1:d], mar=2, distr='zipois', Sigma=Cor, n=n)
    expect_equal(ncol(X), d)
    expect_equal(nrow(X), n)
})


test_that("Synthetic data has correct attributes: negbin", {
    X <- synth_comm_from_counts(amgut1.filt.cs[,1:d], mar=2, distr='negbin', Sigma=Cor, n=n)
    expect_equal(ncol(X), d)
    expect_equal(nrow(X), n)
})

test_that("Synthetic data has correct attributes: pois", {
    X <- synth_comm_from_counts(amgut1.filt.cs[,1:d], mar=2, distr='pois', Sigma=Cor, n=n)
    expect_equal(ncol(X), d)
    expect_equal(nrow(X), n)
})

test_that("Synthetic data has correct attributes: lognorm", {
    X <- synth_comm_from_counts(amgut1.filt.cs[,1:d], mar=2, distr='pois', Sigma=Cor, n=n)
    expect_equal(ncol(X), d)
    expect_equal(nrow(X), n)
})
