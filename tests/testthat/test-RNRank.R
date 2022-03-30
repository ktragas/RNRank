test_that("NULL matrix returns NULL if not throwing", {
  expect_equal(RNRank(NULL,throwOnError = F), NULL)
})
