#' Returns internal package datasets
#'
#' @param species Which organism to return data for. Everything else except "Human" or "Mouse" returns NULL.
#' @param what `TF` for reference dataset, `miRNA` for microRNA or `GT` for gene table
#'
#' @return Requested dataset or NULL
#' @export
#'
#' @examples
getPackageData<-function(species=c("","Human","Mouse"),what=c("TF","miRNA","GT"))
{
  species=match.arg(species)
  if (species=="")
    return(NULL)
  what=match.arg(what)
  if (species=="Human") {
    if (what=="TF")
      return(TF_human)
    else if (what=="GT")
      return(GT_human)
    else if (what=="miRNA")
      return(miRNA_human)
    else
      return(NULL)
  }  else { # if (species=="Mouse")
    if (what=="TF")
      return(TF_mouse)
    else if (what=="GT")
      return(GT_mouse)
    else if (what=="miRNA")
      return(miRNA_mouse)
    else
      return(NULL)
  }
}
