## code to prepare `REFERENCE` dataset goes here
TF_human=miRNA_human=kegg_human=keggcat_human=gos_human=
TF_mouse=miRNA_mouse=kegg_mouse=keggcat_mouse=gos_mouse=NULL
TF=miRNA=kegg=keggcat=gos=NULL
RefPath=system.file("extdata","ReferenceFiles/",package="RNRank")
for (species in c("Human","Mouse")) {
  TF_ref=file.path(RefPath,paste0(species,"_TF_plus_oreganno_trrust.tsv"));
  miRNA_ref=file.path(RefPath,paste0(species,"_miRNAs.tsv"));
  kegg_ref=file.path(RefPath,paste0(species,"_kegg.tsv",sep=""));
  keggcat_ref=file.path(RefPath,paste0(species,"_keggcat.tsv"));
  GOs_ref=file.path(RefPath,paste0(species,"_GOs.tsv"));

  TF=utils::read.table(TF_ref,sep="\t",fill=T);
  rownames(TF)=TF[,2];
  miRNA=utils::read.table(miRNA_ref,sep="\t",fill=T);
  rownames(miRNA)=miRNA[,2];
  kegg=utils::read.table(kegg_ref,sep="\t",fill=T);
  rownames(kegg)=kegg[,2];
  keggcat=utils::read.table(keggcat_ref,sep="\t",fill=T);
  rownames(keggcat)=keggcat[,2];
  gos=utils::read.table(GOs_ref,sep="\t",fill=T);
  rownames(gos)=gos[,2];
  if (species=="Human") {
    TF_human=TF
    miRNA_human=miRNA
    kegg_human=kegg
    keggcat_human=keggcat
    gos_human=gos
  } else {  # species=="Mouse"
    TF_mouse=TF
    miRNA_mouse=miRNA
    kegg_mouse=kegg
    keggcat_mouse=keggcat
    gos_mouse=gos
  }
} # for (species...
usethis::use_data(TF_human, miRNA_human,kegg_human,keggcat_human,gos_human,
                  TF_mouse, miRNA_mouse,kegg_mouse,keggcat_mouse,gos_mouse,
                  internal=T, overwrite = T, compress = "xz")
