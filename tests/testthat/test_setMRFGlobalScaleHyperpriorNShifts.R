context("Tests setting of MRF hyperpriors")

test_that("MRF settings match known values", {
  expect_equal(0.002093737,
               setMRFGlobalScaleHyperpriorNShifts(100, "HSMRF"))
  expect_equal(0.009376335,
               setMRFGlobalScaleHyperpriorNShifts(100, "GMRF"))
})
