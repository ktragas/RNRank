test_that("RNEAv3 works", {
  expect_output(
    {
      DataPath=system.file("extdata",package="RNRank")
      file=file.path(DataPath,
                     "Input/gene_exp.tsv"
                     #"Input/GSE63889_Diff_0U_120U_gene_exp.diff"
      )
      OutputPath=file.path(rprojroot::find_root(rprojroot::has_dir("tests")),"Output")
      RNEAv3(
        filename=file,
        species="Mouse",network="regulatory",type_of_output="csv",
        output_dir=OutputPath,
        output="gene_exp_out")
      #output="GSE63889")
    }
  )
})
