## code to prepare `REFERENCE` dataset goes here
TF_human=miRNA_human=TF_mouse=miRNA_mouse=NULL
GT_human=GT_mouse=NULL
TF=miRNA=NULL
RefPath=system.file("extdata","ReferenceFiles/",package="RNRank")
for (species in c("Human","Mouse")) {
  # TF_ref=file.path(RefPath,paste0(species,"_TF_plus_oreganno_trrust.tsv"));
  TF_ref=file.path(RefPath,paste0(species,"_TF_Reference.tsv"));
  miRNA_ref=file.path(RefPath,paste0(species,"_miRNAs.tsv"));

  # Biomart gene table
  # Gene stable ID with version, Gene stable ID, Gene Name
  GT_ref=file.path(RefPath,paste0(species,"_GT_mart_export.txt.gz"))

  TF=read.table(TF_ref,sep="\t",fill=T);
  rownames(TF)=TF[,2];
  miRNA=read.table(miRNA_ref,sep="\t",fill=T);
  rownames(miRNA)=miRNA[,2];
  GeneTable=read.table(GT_ref,header=T,sep=",")
  rownames(GeneTable)=GeneTable[,1]
  if (species=="Human") {
    TF_human=TF
    miRNA_human=miRNA
    GT_human=GeneTable
  } else {  # species=="Mouse"
    TF_mouse=TF
    miRNA_mouse=miRNA
    GT_mouse=GeneTable
  }
} # for (species...
usethis::use_data(TF_human, miRNA_human,TF_mouse, miRNA_mouse,
                  GT_human, GT_mouse, internal=T, overwrite = T,
                  compress = "xz")
