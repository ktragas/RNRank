test_that("RNEAv3 works", {
  expect_output(
    {
      file=system.file("extdata",
                       #"Input/gene_exp.tsv",
                       "Input/GSE63889_Diff_0U_120U_gene_exp.diff",
                       package="RNRank")
      setwd(system.file(package = "RNRank"))

      RNEAv3(
        filename=file,
        species="Mouse",network="regulatory",type_of_output="csv",
        output_dir="Output",
        #output="gene_exp_out")
      output="GSE63889")
    }
  )
})

test_that("NULL matrix", {
  expect_equal(RNRank(NULL,throwOnError = F), NULL)
  expect_error(RNRank(NULL))
})
