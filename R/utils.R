isValidNetworkMatrix<-function(srcm, throwOnError=F) {
  isValid=(is.matrix(srcm) && ncol(srcm)>=2 && nrow(srcm)>=2 && typeof(srcm)=="character")
  if (!isValid && throwOnError)
    stop("Need character matrix with at least 2 rows and 2 columns. Can't continue")
  return(isValid)
}
