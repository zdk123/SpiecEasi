baseline <- readRDS(test_path("fixtures", "slr_baseline.rds"))

describe("SLR baseline (PR #292 solver)", {
  for (nm in names(baseline$cases)) {
    test_that(paste("fit reproduces the baseline fixture:", nm), {
      fit <- slr_fit_case(baseline$cases[[nm]])
      ref <- baseline$fits[[nm]]
      expect_equal(fit$lambda, ref$lambda)
      for (i in seq_along(ref$icov)) {
        expect_lt(rel_err(fit$icov[[i]],  ref$icov[[i]]),  1e-8)
        expect_lt(rel_err(fit$resid[[i]], ref$resid[[i]]), 1e-8)
      }
    })
  }
})

describe("Partial eigen L-update (dsyevr)", {
  sym <- function(p, seed) { set.seed(seed); X <- matrix(rnorm(p * p), p); (X + t(X)) / 2 }
  spec <- function(d, seed) {
    set.seed(seed); Q <- qr.Q(qr(matrix(rnorm(length(d)^2), length(d))))
    Q %*% diag(d) %*% t(Q)
  }

  test_that("matches the full eigendecomposition on indefinite inputs", {
    for (p in c(10, 50, 200)) {
      M <- sym(p, p)
      lmax <- max(eigen(M, symmetric=TRUE, only.values=TRUE)$values)
      for (k in unique(c(0, 1, 3, p - 1, p)))
        expect_lt(rel_err(.eigPSDThreshPartial(M, k, 0), .eigPSDThresh2(M, k, 0)), 1e-10)
      for (beta in c(.1, 1, lmax + 1))
        expect_lt(rel_err(.eigPSDThreshPartial(M, 2, beta), .eigPSDThresh2(M, 2, beta)), 1e-10)
    }
  })

  test_that("is exact on tied and near-tied spectra", {
    d <- c(5, 5, 5, 2, seq(1, -3, length.out=36))       # three-way tie in the top k
    M <- spec(d, 1)
    for (k in 1:4) expect_lt(rel_err(.eigPSDThreshPartial(M, k, 0), .eigPSDThresh2(M, k, 0)), 1e-10)
    d2 <- c(4, 3, 2 + 1e-8, 2, seq(1, -3, length.out=36)) # near-tie at the k/k+1 cutoff
    M2 <- spec(d2, 2)
    expect_lt(rel_err(.eigPSDThreshPartial(M2, 3, 0), .eigPSDThresh2(M2, 3, 0)), 1e-10)
  })

  test_that("handles eigenvalues spread over many orders of magnitude", {
    d <- c(10^seq(4, -8, length.out=30), -10^seq(-8, 2, length.out=10))
    M <- spec(d, 3)
    for (k in c(1, 3, 10)) {
      err <- norm(.eigPSDThreshPartial(M, k, 0) - .eigPSDThresh2(M, k, 0), 'F') / max(abs(d))
      expect_lt(err, 1e-12)
    }
  })

  for (nm in names(baseline$cases)) {
    test_that(paste("full fit with partial.eig reproduces the baseline:", nm), {
      fit <- slr_fit_case(baseline$cases[[nm]], opts=list(partial.eig=TRUE))
      ref <- baseline$fits[[nm]]
      for (i in seq_along(ref$icov)) {
        expect_lt(rel_err(fit$icov[[i]],  ref$icov[[i]]),  1e-8)
        expect_lt(rel_err(fit$resid[[i]], ref$resid[[i]]), 1e-8)
      }
    })
  }
})

describe("R-update from one eigendecomposition", {
  # Analytic reference: K = Q diag(d) Q' with known Q and d, and the prox
  # (sqrt(K^2 + 4 mu I) - K)/2 evaluated per eigenvalue without cancellation.
  prox_ref <- function(d, Q, mu) {
    s <- sqrt(d^2 + 4 * mu)
    t <- ifelse(d >= 0, 2 * mu / (s + d), (s - d) / 2)
    Q %*% diag(t) %*% t(Q)
  }
  basis <- function(p, seed) { set.seed(seed); qr.Q(qr(matrix(rnorm(p * p), p))) }

  spectra <- list(
    well = seq(-3, 3, length.out=40),
    kappa_1e4 = c(10^seq(2, -2, length.out=20), -10^seq(-2, 2, length.out=20)),
    kappa_1e8 = c(10^seq(4, -4, length.out=20), -10^seq(-4, 4, length.out=20)),
    large_d   = c(seq(80, 120, length.out=20), -seq(80, 120, length.out=20))
  )
  for (nm in names(spectra)) for (mu in c(1, 1e-2, 1e-4)) {
    test_that(sprintf("matches the analytic prox: %s, mu = %g", nm, mu), {
      d <- spectra[[nm]]; Q <- basis(length(d), 7)
      K <- Q %*% diag(d) %*% t(Q); K <- (K + t(K)) / 2
      ref <- prox_ref(d, Q, mu)
      new <- .Rupdate_new(K, mu)
      old <- .Rupdate_old(K, mu)
      err_new <- rel_err(new, ref); err_old <- rel_err(old, ref)
      expect_lt(err_new, 1e-12)
      expect_lte(err_new, max(err_old, 1e-14))
      expect_gt(min(eigen(new, symmetric=TRUE, only.values=TRUE)$values), 0)
    })
  }

  for (nm in names(baseline$cases)) {
    test_that(paste("full fit with fast.R agrees with the baseline:", nm), {
      fit <- slr_fit_case(baseline$cases[[nm]], opts=list(fast.R=TRUE))
      ref <- baseline$fits[[nm]]
      for (i in seq_along(ref$icov)) {
        expect_lt(rel_err(fit$icov[[i]],  ref$icov[[i]]),  1e-6)
        expect_lt(rel_err(fit$resid[[i]], ref$resid[[i]]), 1e-6)
      }
    })
  }
})

describe("Path-changing options converge to the same optimum", {
  # These options change the iterates, so they are checked against a tightly
  # converged reference rather than against the baseline at the default tol.
  cs  <- baseline$cases$rank_np
  ref <- slr_fit_case(cs, opts=list(tol=1e-11, maxiter=1e6))
  configs <- list(balance_mu   = list(balance.mu=TRUE),
                  sparse_first = list(path.order='sparse_first'),
                  both         = list(balance.mu=TRUE, path.order='sparse_first'))
  for (nm in names(configs)) {
    test_that(paste("reaches the reference solution:", nm), {
      fit <- slr_fit_case(cs, opts=c(configs[[nm]], list(tol=1e-8, maxiter=1e5)))
      err <- max(mapply(rel_err, fit$icov, ref$icov), mapply(rel_err, fit$resid, ref$resid))
      expect_lt(err, 1e-3)
    })
  }
})

describe("SLR low-rank component", {
  test_that("L is symmetric PSD with rank at most r", {
    fit <- baseline$fits$rank_np
    r   <- baseline$cases$rank_np$r
    for (L in fit$resid) {
      expect_lt(max(abs(L - t(L))), 1e-10)
      ev <- eigen(L, symmetric=TRUE, only.values=TRUE)$values
      expect_gt(min(ev), -1e-8 * max(1, max(ev)))
      expect_lte(sum(ev > 1e-8 * max(1, max(ev))), r)
    }
  })
})
