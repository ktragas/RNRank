test_that("NULL matrix", {
  expect_equal(RNRank(NULL,throwOnError = F), NULL)
  expect_error(RNRank(NULL))
})
