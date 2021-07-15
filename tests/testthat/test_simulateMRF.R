context("Tests simulation of MRF")

test_that("MRF simulation matches", {
  set.seed(42)
  x <-
    simulateMRF(
      n_episodes = 100,
      model = "HSMRF",
      global_scale_hyperprior = 0.0021
    )
  expect_equal(1.0351365237, mean(x))
  expect_equal(1.0628421716, median(x))
  expect_equal(1.0, x[1])
  expect_equal(0.8160889846, x[100])
  expect_equal(0.9739524255, sum(abs(x[-1] - x[-100])))
  expect_equal(0.3520475309, max(x) - min(x))
  expect_equal(99, which.min(x))
  expect_equal(42, which.max(x))
})
