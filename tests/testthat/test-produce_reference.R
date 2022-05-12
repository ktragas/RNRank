test_that("Create reference files", {
  expect_output({
    organisms=data.frame(id=c("10090","9606"),
                         tf_prefix=c("Mus_musculus","Homo_sapiens"),
                         organism=c("Mouse","Human"))

    DataRoot=system.file("extdata",package="RNRank")
    DataPath=file.path(DataRoot,"ReferenceSources")
    OutputPath=file.path(DataRoot,"ReferenceFiles")

    apply(organisms,1,function(x) {
      pi_file=paste0(x["id"],".","protein.links.full.v11.5.txt.gz")
      tf_file=paste0(x["tf_prefix"],"_TF.txt")
      bm_file=paste0(x["organism"],"_mart_export.txt.gz") # BioMart
      # bm_file=paste0(x["id"],".","protein.info.v11.5.txt.gz") # STRING
      pi_file=file.path(DataPath,pi_file)
      tf_file=file.path(DataPath,tf_file)
      bm_file=file.path(DataPath,bm_file)
      out_file=paste0(x["organism"],"_TF_Reference.tsv")
      out_file=file.path(OutputPath,out_file)
      produce_reference_from_files(pi_file,tf_file,2,bm_file,out_file,
              seps=c("","\t",","), pg_format="Normal",
              # headers=c(T,T,F), seps=c("","\t","\t"), pg_format="STRING",
              org_id=x["id"],
              thresholds=c(0,0), thresholdMethod="or",
              as.data.frame=T)
    })
  })
})
