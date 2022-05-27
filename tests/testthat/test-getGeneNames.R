test_that("get known gene names", {
  expect_equal({
      src=c("ENSG00000163788.9","ENSG00000129657.10","ENSG00000186105.7")
      getGeneNames("Human",src)
    },
    {
      c("SNRK","SEC14L1","LRRC70")
    }, ignore_attr=T)
})
