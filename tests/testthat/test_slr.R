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
