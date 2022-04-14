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
