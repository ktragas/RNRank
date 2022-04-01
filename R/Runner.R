# debugSource("RNEA_v2.R")
# debugSource("RNEA_v3.R")
# sink()
# RNEAv3(
#    filename="Input/gene_exp.tsv",
#    #filename="Input/GSE63889_Diff_0U_120U_gene_exp.diff",
#    species="Mouse",network="regulatory",type_of_output="csv",
#    output_dir="Output",
#    output="gene_exp_out")
#    #output="GSE63889")

#debugSource("R/TFRank.R")
#sink()
srcm=as.matrix(read.table(
   "Output/gene_exp_out_Network.csv",
   #"Output/GSE63889_Network.csv",
   header=T,sep=","))
 P=RNRank(srcm, max_iterations = 200, threshold=0.01, damping=0.85,
          self=T, letZeros = T, divider = 1000.0, verbose = T)
 print(head(P,10))
