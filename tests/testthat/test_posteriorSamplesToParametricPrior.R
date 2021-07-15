context("Tests simulation of MRF")

test_that("MRF simulation matches", {
  set.seed(47)
  gp <-
    posteriorSamplesToParametricPrior(rgamma(1e4, 1, 2), "gamma", 2.0)
  np <-
    posteriorSamplesToParametricPrior(rnorm(1e4, -2, 3), "normal", 2.0)

  gp.ref <- c(0.4771022693, 0.9396894192)
  names(gp.ref) <- c("gamma.shape", "gamma.rate")
  np.ref <- c(-1.9965172801, 4.1974666979)
  names(np.ref) <- c("normal.mean", "normal.sd")

  expect_equal(gp.ref[1], gp[1])
  expect_equal(gp.ref[2], gp[2])

  expect_equal(np.ref[1], np[1])
  expect_equal(np.ref[2], np[2])
})
