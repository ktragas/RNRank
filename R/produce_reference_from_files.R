

#' @title  Produce a list of Transcription Factors interactions
#'
#' @description Produces and saves a file with a list of Transcription Factors interactions,
#' by combining information from provided files#'
#'
#' @param pi_file Protein interactions from STRING (fully detailed, with header)
#' @param tf_file Transcription factors (Tab separated, with header)
#' @param tf_column Column containing gene name of TFs
#' @param bm_file Protein-gene matching from Biomart (N:1, comma separated, with header)
#' @param out_file Results will be saved to this tab-separated file
#' @param headers header parameter for [utils::read.table()] for each of the 3 input files
#' @param seps sep parameter for [utils::read.table()] for each of the 3 input files
#' @param ... Arguments to be passed to [produce_reference()]
#'
#' @return Result of [produce_reference()] or NULL in case of error
#' @export
#'
#' @examples
produce_reference_from_files<-function(pi_file,tf_file,tf_column,bm_file,out_file,
                                       headers=rep(T,3), seps=rep("",3),
                                       ...)
{
  headers=as.logical(headers)
  if (!is.logical(headers)) headers=rep(T,3)
  else if ((lh=length(headers))<3) headers=c(headers,rep(T,3-lh))

  seps=as.character(seps)
  if (!is.character(seps)) seps=rep("",3)
  else if ((ls=length(seps))<3) seps=c(seps,rep("",3-ls))

  # Φόρτωση πρωτεϊνικών αλληλεπιδράσεων
  t=read.table(pi_file,header = headers[1],sep=seps[1])

  # Φόρτωση μεταγραφικών παραγόντων
  tf=read.table(tf_file,header=headers[2],sep=seps[2])

  # Φόρτωση συμπιεσμένου CSV αρχείου biomart
  # (περιέχει εγγραφές "Protein stable ID,Gene name")
  biomart=read.table(bm_file, header=headers[3], sep=seps[3])

  ref=produce_reference(t,tf[,tf_column],biomart,...)

  # Σώσιμο αποτελεσμάτων
  if (!is.null(ref)) {
    if (is.data.frame(ref)) {
      write.table(ref,out_file,append = F,sep="\t",quote=F,na="",
                  row.names = F,col.names = F)
    } else { # list
      proclist=function(x) {
        rec=c(x$len,x$tf,x$genes)
        paste0(rec,collapse="\t")
      }
      ref=lapply(ref,proclist)
      unlink(out_file)
      lapply(ref,write,file=out_file,append=T, sep="\t")
    }
  }
  return(ref)
}
