#' Creates full regulatory network from reference files
#'
#' @param x Use only these Transcription Factors. `NULL` (the default) creates
#' the biggest possible network.
#' @param species Organism. Respective reference file must exist.
#' @param suffix Reference file name is constructed by pasting `species` and `suffix`.
#' @param reference_dir Directory containing the reference file. `?` (the default) uses
#' the reference directory of the package.
#' @param verbose Shows some information during process.
#'
#' @return A 2-column matrix representing directed edges of the network
#' @export
#'
#' @examples
makeNetFromReference<-function(x=NULL,species="Mouse",suffix="_TF_Reference.tsv",reference_dir="?",verbose=T)
{
  if (reference_dir=="?") {
    reference_dir=system.file("extdata","ReferenceFiles",package="RNRank")
  }
  TF_ref=file.path(reference_dir,paste0(species,suffix));
  TF=read.table(TF_ref,sep="\t",fill=T);
  rownames(TF)=TF[,2];

  if (!is.null(x)) {
    x=as.vector(x)
    x=intersect(x,rownames(TF))
    TF=TF[x,]
  }

  Network=matrix(ncol=2,nrow=0)
  colnames(Network)=c("Source","Target");
  a=apply(TF,1,function(x) {
    last_index=2+as.integer(x[1])
    return(cbind(x[2],x[3:last_index]))
  })
  Network=do.call("rbind", a)
  rownames(Network)=NULL
  if (verbose) {
    ng=length(unique(c(Network)))
    nr=nrow(Network)
    if (nr>0)
      msg=sprintf("Network consists of %d genes and contains %d edges",ng,nr)
    else
      msg="Empty network"
    print(msg)
  }
  return(Network)
}
