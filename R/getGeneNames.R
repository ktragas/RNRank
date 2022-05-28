#' Replaces Ensembl IDs with gene names
#'
#' @description Tries to replace the Ensembl IDs ("ENSGxxxxx.x") with gene names, using
#' internally stored data.
#'
#' @param species Organism. "Human","Mouse" datasets exist in package.
#' @param x The vector of Ensembl IDs to be replaced with the corresponding
#'
#' @return An equally-sized vector with replaced values, where possible.
#' @export
#'
#' @examples
getGeneNames<-function(species=c("","Human","Mouse"),x)
{
  species=match.arg(species)
  if (species=="")
    return(x)
  else if (species=="Human")
    tbl=GT_human
  else # if (species=="Mouse")
    tbl=GT_mouse
  genes=sapply(x, function(x,species,bm) {
    # Gene stable ID with version, Gene stable ID, Gene Name
    # Ψάχνω να βρω το ακριβές, μετά χωρίς τη version και τέλος ως έχει
    g=bm[x==bm[,1],3]
    if (length(g)==0) {
      g=bm[startsWith(x,bm[,2]),3]
    }
    if (is.na(g) || length(g)==0 || trimws(g)=="") {
      g=x
    }
    return(g)
  }, species,tbl)
  return(genes)
}
