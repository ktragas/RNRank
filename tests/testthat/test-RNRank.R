test_that("NULL matrix", {
  expect_equal(RNRank(NULL,throwOnError = F), NULL)
  expect_error(RNRank(NULL))
})

test_that("RNRank works", {
  expect_output({
    OutputPath=file.path(rprojroot::find_root(rprojroot::has_dir("tests")),"Output")
    RNEA_output_file=file.path(OutputPath,"gene_exp_out_Network.csv")
    srcm=as.matrix(read.table(RNEA_output_file,header=T,sep=","))
    P=RNRank(srcm, max_iterations = 200, threshold=0.01, damping=0.85,
             self=T, letZeros = T, divider = 1000.0, verbose = T)
    print(head(P,10))
    barplot(P[1:min(length(P),30),1],las=2,main="Most important genes")
    rm(RNEA_output_file,srcm,P)
  })
})

