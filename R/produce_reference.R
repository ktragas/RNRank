#' @title  Produce a list of Transcription Factors interactions
#'
#' @description Produces and returns a list of Transcription Factors interactions,
#' by combining provided information.
#'
#' @details Two methods are supported for threshold checking. The default is `and`
#' where the thresholds are checked one after another, and the other is `or`,
#' where each condition is checked separately and the results are combined.
#'
#' @param t Protein interactions from STRING database.
#' @param tf Transcription factors (gene names).
#' @param biomart Protein-Gene matching table. Normally from BioMart database.
#' @param org_id STRING organism ID. If specified, it is assumed that protein names
#' in the interactions table will have the form *<org_id>.<protein>*.
#' @param thresholdCols Column numbers in protein interactions dataset,
#' to be checked against thresholds
#' @param thresholds Respective thresholds for thresholdCols.
#' Only greater values than the threshold are accepted.
#' @param thresholdMethod The method to be used for threshold checking. See 'Details'.
#' @param as.data.frame Return data frame instead of list
#' @param verbose Print stats during processing
#'
#' @return List or data frame containing the TF interactions. NULL in case of error.
#' Format is: Count of interactions, TF, Interactors
#' @export
#'
#' @examples
produce_reference<-function(t,tf,biomart,org_id="",thresholdCols=c(8,10),thresholds=c(),
                            thresholdMethod=c("and","or"),as.data.frame=F,verbose=T)
{
  # Έλεγχος παραμέτρων
  nc=ncol(t)
  if (!is.numeric(nc) || nc<2) return(NULL)
  pi_len=nrow(t)
  if (verbose){
    print(sprintf("%d total protein interactions",pi_len))
  }

  tf=as.vector(tf)
  if (!is.vector(tf)) return(NULL)
  # Διαπίστωσα ότι στο ποντίκι είχε διπλές εγγραφές
  tf=unique(trimws(tf))
  tf=tf[nchar(tf)>0]
  tf_len=length(tf)
  if (verbose) {
     print(sprintf("%d TFs",tf_len))
  }
  if (tf_len<=0) return(NULL)

  org_id=as.character(org_id)
  if (!is.character(org_id) || length(org_id)!=1)
    return(NULL)
  if (nchar(org_id)>0) org_id=paste0(org_id,".")

  lc=length(thresholdCols)
  lt=length(thresholds)
  if (lt<lc) thresholds=c(thresholds,rep(-Inf,lc-lt)) # Αν τα κατώφλια είναι λιγότερα, συμπληρώνω -Inf
  if (lc>0) {
    # Παραλείπω τους μη έγκυρους αριθμούς στήλης και τα αντίστοιχα κατώφλια
    # Αν επίσης κάποιο κατώφλι δεν είναι αριθμητικό γίνεται -Inf
    tmpc=c()
    tmpt=c()
    for (i in 1:lc) {
      if (thresholdCols[i]<=2 || thresholdCols[i]>nc) next  # Οι πρώτες 2 στήλες είναι πάντα οι πρωτεΐνες
      tmpc=c(tmpc,thresholdCols[i])
      tmp=as.numeric(thresholds[i])
      tmpt=c(tmpt,ifelse(is.numeric(tmp),tmp,-Inf))
    }
    thresholdCols=tmpc
    thresholds=tmpt
  }

  # Αφαιρώ γραμμές χωρίς protein id και γραμμές χωρίς gene name
  bm_orglen=nrow(biomart) # Αρχικός αριθμός εγγραφών
  w=which(biomart[,1]=="" | biomart[,2]=="")
  if (length(w)>0)
    biomart=biomart[-w,]
  bm_len=nrow(biomart) # Τελικός αριθμός εγγραφών
  if (verbose) {
    print(sprintf("Keeping %d (out of %d) protein-gene matching records",bm_len,bm_orglen))
  }

  # Μετατρέπω σε named vectors, με κωδικό οργανισμού, για πιο εύκολες αναζητήσεις
  # Προϋπόθεση να μην αντιστοιχεί μία πρωτεΐνη σε πάνω από ένα γονίδιο.
  bm_id=biomart[,2]
  names(bm_id)=paste0(org_id,biomart[,1])

  # Εγγραφές της biomart που τα γονίδιά τους αντιστοιχούν στους μεταγραφικούς παράγοντες
  # Δεν είναι 100% απαραίτητο, το κάνω για να επιταχύνω τα ψαξίματα
  w=which(biomart[,2] %in% tf)
  tf_bm=biomart[w,]

  # Από τις αλληλεπιδράσεις, κρατάω όσες περιέχουν τουλάχιστον μία πρωτεΐνη μεταγραφικού παράγοντα.
  # Καθώς μπορεί να περιέχουν κωδικό οργανισμού, φροντίζω να κάνω κατάλληλο έλεγχο
  if (nchar(org_id)==0) checkCol=1
  else {
    tf_bm[,3]=paste0(org_id,tf_bm[,1])
    checkCol=3
  }
  t=t[(t[,1] %in% tf_bm[,checkCol] | t[,2] %in% tf_bm[,checkCol]),]
  pi_len=nrow(t)
  if (verbose){
    print(sprintf("%d TF protein interactions before filtering",pi_len))
  }

  # Για να λιγοστέψουν τα δεδομένα στη μνήμη, κρατώ μόνο τις στήλες που με ενδιαφέρουν,
  # δηλ. τις πρωτεΐνες και τις thresholdCols
  t1=data.frame(protein1=t[,1],protein2=t[,2])
  for (c in thresholdCols) {
    name=names(t[c])[1]
    t1[name]=t[c]
  }
  t=t1

  # Έλεγχος κατωφλίων
  thresholdMethod=match.arg(thresholdMethod)
  for (i in seq_along(thresholdCols)) {
    # Οι στήλες πια είναι οι 2 πρωτεΐνες και όσες thresholdCols είχαν επιλεγεί,
    # με τη σειρά που είχαν επιλεγεί
    tsel=t[t[,2+i]>thresholds[i],]
    if (thresholdMethod=="or") {
        if (i==1) {
          lc=length(thresholdCols)
          t1=NULL
        }
        t1=rbind(t1,tsel)
        if (i==lc)
          t=unique(t1)
    }
    else   # Default: thresholdMethod=="and"
      t=tsel
  }
  rm(t1)
  pi_len=nrow(t)
  if (verbose){
    print(sprintf("%d TF protein interactions after filtering",pi_len))
  }
  # Για κάθε μεταγραφικό παράγοντα, βρίσκω τις συνδέσεις του και φτιάχνω μια λίστα
  # που αντιστοιχεί σε μία εγγραφή του αρχείου Reference. Η λίστα ref κρατάει όλες
  # αυτές τις λίστες.
  ref=list()
  li=1
  for (i in seq_along(tf)) {
    print(sprintf("[%d] %s",i,tf[i]))
    # Ποιες πρωτεΐνες έχει η biomart για αυτόν τον παράγοντα (γονίδιο)
    w=which(tf_bm[,2]==tf[i])
    tf_proteins=tf_bm[w,checkCol]

    # Ποιες πρωτεΐνες αλληλεπιδρούν με αυτές
    cons=unique(c(t[which(t[,1] %in% tf_proteins),2],
                  t[which(t[,2] %in% tf_proteins),1]))
    if (length(cons)==0) next

    cons=unique(bm_id[cons]) # Μετατροπή σε γονίδια
    # Κάποιες πρωτεΐνες δεν αντιστοιχούν σε γνωστό γονίδιο και αφαιρούνται
    cons=cons[!(is.na(cons) | is.null(cons))]
    len=length(cons)
    if (len==0) next

    # Έχω όλα τα στοιχεία του αρχείου Reference
    record=list(count=length(cons), tf=tf[i], genes=sort(cons))
    ref[[li]]=record
    li=li+1
  }

  # Ταξινόμηση κατά πλήθος ρυθμιζόμενων
  lengths=sapply(ref,"[[",1) # ή sapply(ref,"[[","count")
  ref=ref[order(lengths,decreasing = T)]

  if (verbose) {
     print(sprintf("%d Reference records with %d maximum single TF interactions",
                 length(ref), ref[[1]]$count))
  }

  # Μετατροπή σε data.frame
  if (as.data.frame==T) {
    ref=do.call(dplyr::bind_rows,lapply(ref,function(x) {
      # Όταν το x$genes είναι μόνο ένα, δεν παίρνει σωστό όνομα η στήλη και προστίθεται
      # πολύ μακριά με την bind_rows(). Οπότε προσθέτω ένα κενό μέλος
      if (length(x$genes)==1) x$genes=c(x$genes,"")
      data.frame(x$count,x$tf,t(x$genes))
    }))
    ref[is.na(ref)]=""
  }
  return(ref)
}
