## Regenerate the SLR baseline fixture used by test_slr.R.
## The stored fixture was produced by the PR #292 solver (commit 1dc723d on perf/slr-admm).
## Bit-compatible solver changes must reproduce it; run from the package root:
##   Rscript tests/testthat/fixtures/make_slr_baseline.R

pkgload::load_all(".", quiet=TRUE)
source("tests/testthat/helper-slr.R")

cases <- slr_baseline_cases()
fits  <- lapply(cases, function(cs) slr_fit_case(cs))
saveRDS(list(cases=cases, fits=fits), "tests/testthat/fixtures/slr_baseline.rds", version=2)
