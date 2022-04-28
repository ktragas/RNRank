test_that("RNEAv3 works", {
  expect_output(
    {
      DataPath=system.file("extdata",package="RNRank")
      file=file.path(DataPath,
                     #"Input/gene_exp.tsv"
                     "Input/GSE63889_Diff_0U_120U_gene_exp.diff"
                     )
      ReferenceDir=file.path(DataPath,"ReferenceFiles")
      OutputPath=file.path(rprojroot::find_root(rprojroot::has_dir("tests")),"Output")
      RNEAv3(
        filename=file,species="Mouse",internal_data=F,
        network="regulatory",type_of_output="csv",
        reference_dir=ReferenceDir,output_dir=OutputPath,
        #output="gene_exp_out"
        output="GSE63889"
        ,rank=T,max_iterations = 200, threshold=0.001, damping=0.85,
        self=T, letZeros = T, divider = 1000.0, verbose = T
        )

      #unlink(OutputPath,recursive = T,force = T)
    }
  )
})
