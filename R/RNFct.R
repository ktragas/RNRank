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
  # ----------------------------
  return(RNFct(network,"ffl",verbose,throwOnError))
} # end function

#' Locate circles in regulatory network
#'
#' @description Locate circles in regulatory network.
#'
#' @eval net_param()
#' @eval common_params()
#'
#' @return 3-column matrix of circles. Column X regulates column Y, Y regulates Z
#' and Z regulates X.
#' @export
#'
#' @examples
RNCircle<-function(network,verbose=F,throwOnError=T)
{
  # Εντοπισμός κύκλου στο υποδίκτυο
  #        X---->Y---->Z
  #         \----<----/
  # ------------------------------
  return(RNFct(network,"circle",verbose,throwOnError))
} # end function


RNFct<-function(network,func=c("ffl","circle"),verbose=F,throwOnError=T)
{
  # Έλεγχος παραμέτρων
  if (!isValidNetworkMatrix(network,throwOnError))
    return(NULL)

  # Αφαιρώ διπλές εγγραφές και εγγραφές αυτορρύθμισης
  # και αρχικοποιώ πίνακα FCT
  network=unique(network)
  ngenes=length(unique(c(network)))
  sn=which(network[,1]==network[,2])
  if (length(sn)>0)
    network=network[-sn,]
  fct=matrix(nrow=0,ncol=3)
  colnames(fct)=c("X","Y","Z")

  func=match.arg(func)

  if (func=="ffl") {
    title="FFL"

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
        fct=rbind(fct,cbind(x,y,z))
      } # for (z...
    } # for (x...
  }
  else # if (func=="circle")
  {
    title="circle"

    # Σε κύκλο μπορούν να συμμετέχουν μόνο γονίδια που είναι και ρυθμίζοντα
    # και ρυθμιζόμενα, άρα οι ακμές που περιέχουν τα υπόλοιπα αφαιρούνται
    cm=intersect(network[,1],network[,2])
    network=network[network[,1] %in% cm & network[,2] %in% cm, ]
    # Εντοπίζω όλα τα γονίδια από τα οποία ξεκινάει ακμή και τα θεωρώ
    # αρχή (x) πιθανού κύκλου
    for (x in unique(network[,1])) {
      # Τα γονίδια που ρυθμίζονται από το x είναι πιθανά y
      # και αυτά που ρυθμίζουν το x είναι πιθανά z
      possibleY=unique(network[network[,1]==x,2])
      possibleZ=unique(network[network[,2]==x,1])

      for (z in possibleZ) {
        # Από τα ρυθμιζόμενα από το x, όσα ρυθμίζουν το z, είναι y
        y=network[network[,2]==z & network[,1] %in% possibleY, 1]
        if (length(y)<=0) next

        # Βρέθηκαν!
        fct=rbind(fct,cbind(x,y,z))
      } # for (z...
    } # for (x...
  }
  rownames(fct)=NULL
  if (verbose) {
    n=nrow(fct)
    if (n>0)
      msg=sprintf("Found %d %ss consisting of %d (out of %d) different genes",
                  n,title,length(unique(c(fct))),ngenes)
    else
      msg=sprintf("No %ss found",title)
    print(msg)
  }
  return(fct)
}
