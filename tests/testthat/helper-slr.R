## Shared setup for the SLR tests and the baseline fixture.

slr_sim_data <- function(p, n, r=2, seed=10010) {
  set.seed(seed)
  graph <- make_graph('band', p, p)
  Prec  <- graph2prec(graph)
  Cor   <- cov2cor(prec2cov(Prec))
  Z     <- rmvnorm(n, rep(0, p), Cor)
  U     <- matrix(rnorm(n * r), n, r)
  W     <- matrix(rnorm(p * r, sd=.5), p, r)
  Z + U %*% t(W)
}

slr_baseline_cases <- function() {
  list(
    rank_np   = list(p=20, n=60, r=2, beta=NULL, nlambda=5),
    rank_pn   = list(p=30, n=15, r=2, beta=NULL, nlambda=5),
    penalty   = list(p=20, n=60, r=NULL, beta=.1, nlambda=5)
  )
}

## Fixture cases carry their data: rmvnorm goes through eigen/svd, so the simulated X
## differs across LAPACK builds (eigenvector signs).
slr_fit_case <- function(cs, opts=NULL) {
  X <- if (!is.null(cs$X)) cs$X else slr_sim_data(cs$p, cs$n)
  args <- list(data=X, cor=TRUE, nlambda=cs$nlambda, lambda.min.ratio=1e-2,
               lambda.max=getMaxCov(cor(X)), opts=opts)
  if (!is.null(cs$r))    args$r    <- cs$r
  if (!is.null(cs$beta)) args$beta <- cs$beta
  fit <- do.call(sparseLowRankiCov, args)
  list(icov  = lapply(fit$icov,  as.matrix),
       resid = lapply(fit$resid, as.matrix),
       lambda = fit$lambda)
}

rel_err <- function(a, b) norm(as.matrix(a) - as.matrix(b), 'F') / max(1, norm(as.matrix(b), 'F'))
