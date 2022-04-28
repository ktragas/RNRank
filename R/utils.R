isValidNetworkMatrix<-function(srcm, throwOnError=F) {
  isValid=(is.matrix(srcm) && ncol(srcm)>=2 && nrow(srcm)>=2 && typeof(srcm)=="character")
  if (!isValid && throwOnError)
    stop("Need character matrix with at least 2 rows and 2 columns. Can't continue")
  return(isValid)
}

net_param<-function() {
  c("@param network Two-column gene matrix. Every row contains one regulating and one regulated gene.")
}

common_params<-function() {
  c(
    "@param verbose Shows some information during process.",
    "@param throwOnError If FALSE, in case of stopping error, returns NULL. Otherwise [base::stop()] is called."
  )
}

makeNetFromReference<-function(species="Mouse",reference_dir="?")
{
   if (reference_dir=="?") {
     reference_dir=system.file("extdata","ReferenceFiles",package="RNRank")
   }
   TF_ref=file.path(reference_dir,paste0(species,"_TF_Reference.tsv"));
   TF=read.table(TF_ref,sep="\t",fill=T);
   rownames(TF)=TF[,2];
   Network=matrix(ncol=2,nrow=0)
   colnames(Network)=c("Source","Target");
   Network=do.call("rbind", apply(TF,1,function(x) {
     last_index=2+as.integer(x[1])
     rbind(Network,cbind(x[2],x[3:last_index]))
   }))
   rownames(Network)=NULL
   return(Network)
}
