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

context("marginal distribuitions applied correctly")
# This test has been structured to mimic the relevant functionality,
# since currently we do not have a way to access the normal samples used internally.
test_that("VGAM::qzinegbin is applied correctly to normal samples", {
    
    # setup, mimics internal functioning of synth_comm_from_counts and rmvzingbin
    cors <- cor(amgut1.filt.cs)
    distr <- "zinegbin"
    params <- get_comm_params(comm=amgut1.filt.cs, distr=distr)

    paramat <- do.call('rbind', params)
    paramat <- data.frame(apply(paramat, 2, as.numeric))
    
    # Follows the implementation of rmvzingbin
    normal_samps <- rmvnorm(n=n, mu=rep(0, ncol(amgut1.filt.cs)), Sigma=cors)
    unif_probs <- pnorm(normal_samps)

    seRes <- matrix(VGAM::qzinegbin(t(unif_probs), size=paramat$size, 
				    munb=paramat$munb, pstr0=paramat$pstr0), 
		    n, ncol(amgut1.filt.cs), byrow=TRUE)

    col3res <- VGAM::qzinegbin(unif_probs[,3], size=paramat$size[3], 
			       munb=paramat$munb[3], pstr0=paramat$pstr0[3])

    # show(seRes[,3] == col3res)
    expect( all(seRes[,3] == col3res), "After application of the quantile function, did we get the correct counts back?")
})
