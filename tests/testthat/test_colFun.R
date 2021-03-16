context("Tests colFun()")

test_that("colFun() doesn't error out", {
  expect_error(colFun(1), NA)
})
