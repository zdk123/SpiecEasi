## Regenerate the SLR baseline fixture used by test_slr.R.
## The stored fixture was produced by the PR #292 solver with all perf flags off,
## after the double-precision lambda fix in SOFTTHRESH2 (which moved results by ~5e-8).
## Bit-compatible solver changes must reproduce it; run from the package root:
##   Rscript tests/testthat/fixtures/make_slr_baseline.R

pkgload::load_all(".", quiet=TRUE)
source("tests/testthat/helper-slr.R")

cases <- lapply(slr_baseline_cases(), function(cs) c(cs, list(X=slr_sim_data(cs$p, cs$n))))
fits  <- lapply(cases, function(cs) slr_fit_case(cs))
saveRDS(list(cases=cases, fits=fits), "tests/testthat/fixtures/slr_baseline.rds", version=2)
