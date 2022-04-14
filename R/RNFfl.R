#' Locate feed-forward loops in regulatory network
#'
#' @description Locate feed-forward loops in regulatory network.
#'
#' @eval net_param()
#' @eval common_params()
#'
#' @return 3-column matrix of FFLs. Column X regulates columns Y and Z.
#' Column Y regulates Z, too.
#' @export
#'
#' @examples
RNFfl<-function(network,verbose=F, throwOnError=T)
{
  # Εντοπισμός FFL στο υποδίκτυο
  #        X---->Y---->Z
  #         \----->---/
  # ------------------------------

  # Έλεγχος παραμέτρων
  if (!isValidNetworkMatrix(network,throwOnError))
    return(NULL)

  # Αφαιρώ διπλές εγγραφές και εγγραφές αυτορρύθμισης
  # και αρχικοποιώ πίνακα FFL
  network=unique(network)
  ngenes=length(unique(c(network)))
  sn=which(network[,1]==network[,2])
  network=network[-sn,]
  ffl=matrix(nrow=0,ncol=3)
  colnames(ffl)=c("X","Y","Z")

  # Εντοπίζω όλα τα γονίδια από τα οποία ξεκινάει ακμή και τα θεωρώ
  # αρχή (x) πιθανού FFL
  for (x in unique(network[,1])) {
    # Βρίσκω όλα τα ρυθμιζόμενα από το x και τα θεωρώ πιθανό τέλος (z) του FFL
    xtargets=network[network[,1]==x,2]

    # Πρέπει να ρυθμίζει τουλάχιστον 2 γονίδια για να μπορεί να είναι x
    if (length(xtargets)<=1) next

    for (z in xtargets) {
      # Οποιοδήποτε άλλο ρυθμιζόμενο από το x πλην του z, είναι πιθανό y
      possibleY=xtargets[xtargets!=z]

      # Από αυτά, όσα ρυθμίζουν το z, είναι y
      y=network[network[,2]==z & network[,1] %in% possibleY,1]

      # ### Άλλος τρόπος
      # # Όσα άλλα εκτός του x ρυθμίζουν το z, θεωρούνται πιθανά y
      # zparents=network[network[,2]==z,1]
      # zparents=zparents[zparents!=x]
      #
      # # Αυτά που ρυθμίζονται από το x είναι όντως y
      # y=zparents[zparents %in% xtargets]

      if (length(y)<=0) next

      # Βρέθηκαν!
      ffl=rbind(ffl,cbind(x,y,z))
    } # for (z...
  } # for (x...
  rownames(ffl)=NULL
  if (verbose) {
    n=nrow(ffl)
    if (n>0)
      msg=sprintf("Found %d FFLs consisting of %d (out of %d) different genes",
                  n,length(unique(c(ffl))),ngenes)
    else
      msg("FFLs not found")
    print(msg)
  }
  return(ffl)
} # end function
