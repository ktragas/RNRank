test_that("RNEAv3 works", {
  expect_output(
    {
      DataPath=system.file("extdata",package="RNRank")
      file=file.path(DataPath,
                     #"Input/gene_exp.tsv"
                     #"Input/GSE63889_Diff_0U_120U_gene_exp.diff"
                     "Input/GSE182432_RNA-seq_processed_data_file.diff"
                     )
      ReferenceDir=file.path(DataPath,"ReferenceFiles")
      OutputPath=file.path(rprojroot::find_root(rprojroot::has_dir("tests")),"Output")
      output=#"gene_exp_out"
            #"GSE63889"
            "GSE182432"
      #species="Mouse",
      species="Human"
      print(system.time(
        RNEAv3(
        filename=file,species=species, internal_data=F,
        reference_dir=ReferenceDir,output_dir=OutputPath,
        output = output, FC_threshold = 1,ffl = T, circles = T,
        rank=T,max_iterations = 200, threshold=0.001, damping=0.85,
        self=T, letZeros = T, divider = 1000.0, favorings=1.1,
        verbose = T
        )
      ))

      #unlink(OutputPath,recursive = T,force = T)
    }
  )
})
