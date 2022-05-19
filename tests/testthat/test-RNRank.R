test_that("NULL matrix", {
  expect_equal(RNRank(NULL,throwOnError = F), NULL)
  expect_error(RNRank(NULL))
})

# Matrix from http://www.ams.org/publicoutreach/feature-column/fcarc-pagerank
test_that("RNRank with known matrix", {
    expect_equal({
      srcm=matrix(
        data=c("A","B","A","C","B","D","C","B","C","E",
               "D","B","D","E","D","F","E","F","E","G","E","H",
               "F","H","G","A","G","E","G","H","H","F","H","G"),
        ncol=2,byrow = T)
      P=RNRank(srcm, damping=1, letZeros = T, sorted = F)
    },
    {
      m=matrix(data=c(0.06,0.0675,0.03,0.0675,0.0975,0.2025,0.18,0.295),ncol=1)
      colnames(m)=c("Rank")
      rownames(m)=LETTERS[1:8]
      m
    }, tolerance = 0.001 )
})

test_that("RNRank works", {
  expect_output({
    OutputPath=file.path(rprojroot::find_root(rprojroot::has_dir("tests")),"Output")
    name=#"gene_exp_out"
         #"GSE63889"
         #"GSE182432"
        "Human_lung"
    title=sprintf("%s - Most important genes", name)
    RNEA_output_file=file.path(OutputPath,paste0(name,"_Network.csv"))
    srcm=as.matrix(read.table(RNEA_output_file,header=T,sep=","))
    FFL_file=file.path(OutputPath,paste0(name,"_FFLs.csv"))
    ffl=as.matrix(read.table(FFL_file,header=T,sep=","))
    Circle_file=file.path(OutputPath,paste0(name,"_Circles.csv"))
    circles=as.matrix(read.table(Circle_file,header=T,sep=","))
    l=RNRank(srcm, max_iterations = 200, threshold=0.001, damping=0.85,
             self=T, letZeros = T, divider = 1000.0,
             favoredFFLs=ffl,favoredCircles=circles,favorings=c(1.15,1.1),
             returned="both",verbose = T)
    P=l$values
    print(head(P,10))
    # barplot(P[1:min(nrow(P),30),1],las=2,main=title)
    barplot(P[min(nrow(P),20):1,1],las=2,main=title,horiz=T,
            cex.names=0.8,width=10,space=0.4)
  })
})

