context('spiec.easi getters')

load(system.file('testdata', 'sepulsar.RData', package='SpiecEasi'))

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
  expect_equal(getOptInd(out), (i<-out$select$stars$opt.index))
  expect_equal(getOptNet(out), out$refit$stars)
  expect_equal(getRefit(out), out$refit$stars)
  expect_equal(getOptLambda(out), out$lambda[i])
  expect_equal(getOptMerge(out), out$select$stars$merge[[i]])
  expect_equal(getOptBeta(out), out$est$beta[[i]])
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
