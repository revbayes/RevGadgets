context("tests the processSSE function")

test_that("processes SSE traces", {

  # read in and process file
  bisse_file <- system.file("extdata",
                            "sse/primates_BiSSE_activity_period_mini.p",
                            package="RevGadgets")
  pdata <- processSSE(bisse_file)

  # test file format
  expect_equal(class(pdata), "data.frame")

  # check data
  expect_equal(ncol(pdata), 6)
  expect_equal(colnames(pdata), c("value", "rate", "hidden_state", "label",
                                  "observed_state", "Iteration" ))

})
