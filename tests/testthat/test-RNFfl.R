test_that("RNRank works", {
  expect_output({
    OutputPath=file.path(rprojroot::find_root(rprojroot::has_dir("tests")),"Output")
    RNEA_output_file=file.path(OutputPath,
        "GSE63889_Network.csv"
        # "gene_exp_out_Network.csv"
    )
    srcm=as.matrix(read.table(RNEA_output_file,header=T,sep=","))
    P=RNFfl(srcm,verbose = T, throwOnError=F)
    print(head(P,10))
    rm(OutputPath,RNEA_output_file,srcm,P)
  })
})

