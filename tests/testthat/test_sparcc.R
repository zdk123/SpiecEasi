context('setup')

p <- 20
e <- p
n <- 500
set.seed(10010)
g <- make_graph('erdos_renyi', p, e)
S <- cov2cor(prec2cov(graph2prec(g)))
X <- exp(rmvnorm(n, rep(0,p), S))
X.f <- t(apply(X, 1, norm_to_total))

context("Input data")
test_that("data is counts", {
  expect_silent(.data.checks(X))
  expect_warning(.data.checks(X.f))
  expect_warning(.data.checks(scale(X)))
})

context("SparCC")
test_that("sparcc gives expected output", {
  out <- sparcc(X)
  expect_true(all(c('Cov', 'Cor') %in% names(out)))
  expect_true(all(diag(out$Cor)==1))
  expect_false(any(out$Cor==out$Cov))
  expect_equal(out$Cor, cov2cor(out$Cov))
  expect_equal(dim(out$Cor), c(p,p))
  expect_equal(dim(out$Cov), c(p,p))
})

test_that('sparccboot gives expected output', {
  out <- sparccboot(X, R=3)
  expect_equal(dim(out$t), c(3, p*(p-1)/2))
  expect_equal(length(out$t0), p*(p-1)/2)
  expect_equal(dim(out$null_av$t), c(3, p*(p-1)/2))
  expect_equal(length(out$null_av$t0), p*(p-1)/2)

  pvals <- pval.sparccboot(out)
  expect_gte(min(pvals$pvals, na.rm=TRUE), 0)
  expect_lte(max(pvals$pvals, na.rm=TRUE), 1)
  expect_gte(min(pvals$cors,  na.rm=TRUE), -1)
  expect_lte(max(pvals$cors,  na.rm=TRUE), 1)
})
