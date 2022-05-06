test_that("RNFfl works", {
  expect_output({
    name="gene_exp_out"
        #"GSE63889"
        #"GSE182432"
    OutputPath=file.path(rprojroot::find_root(rprojroot::has_dir("tests")),"Output")
    RNEA_output_file=file.path(OutputPath,paste0(name,"_Network.csv"))
    srcm=as.matrix(read.table(RNEA_output_file,header=T,sep=","))
    P=RNFfl(srcm,verbose = T, throwOnError=F)
    print(head(P,10))
    rm(OutputPath,RNEA_output_file,srcm,P,name)
  })
})

