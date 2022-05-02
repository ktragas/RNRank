#' RNRank: A package for gene ranking in inferred active regulatory network
#'
#' The RNRank package provides five  functions:
#' [RNEAv3()], [RNRank()], [RNFfl()], [produce_reference()] and [produce_reference_from_files()]
#'
#' @docType package
#' @keywords internal
#' @importFrom utils read.table write.csv write.table
#' @aliases RNRank-package
"_PACKAGE"

#' Rank by gene importance
#'
#' @description Infers importance of Regulatory Network nodes (genes),
#' by using a Google's PageRank-like algorithm.
#'
#' @eval net_param()
#' @param max_iterations Maximum number of iterations, if not converging earlier (minimum 10).
#' @param threshold Euclidean distance between iterations less than or equal to this value
#' will terminate calculations
#' @param damping Damping factor (0-1) defining percentage of non-randomness.
#' @param self Self regulations permitted (`TRUE`) or not (`FALSE`).
#' @param letZeros Set to `FALSE` to allow random regulations from regulating genes
#' to those not normally regulated.
#' @param divider Divides the percentage used for dangling genes,
#' to use it for random regulations from regulating genes (minimum 10).
#' @param sorted Results sorted by descending percentage (T) or unsorted (F).
#' @param preorder Sort input matrix, so that regulating genes with more targets are first.
#' @eval common_params()
#'
#' @return Named 1-column matrix of percentages (0-1).
#' @export
#'
#' @examples
RNRank = function(network, damping=0.85, max_iterations=100, threshold=0,
                  self=F, letZeros=F, divider=100.0, sorted=T, preorder=F,
                  verbose=F, throwOnError=T)
{
  # Έλεγχος παραμέτρων
  if (!isValidNetworkMatrix(network,throwOnError))
    return(NULL)

  warn=function(x) {
    msg=sprintf("Invalid value for \"%s\". Adjusted to %d",deparse(substitute(x)),x)
    warning(msg)
  }

  if (damping<0 || damping>1) {
    if (damping<0)
      damping=0
    else #if (damping>1)
      damping=1
    warn(damping)
  }
  if (max_iterations<10) {
    max_iterations=100
    warn(max_iterations)
  }
  if (threshold<0) {
    threshold=0
    warn(threshold)
  }
  if (divider<10) {
    divider=100
    warn(divider)
  }

  # Αφαιρώ γονίδια που αυτορρυθμίζονται (source==target)
  m=network
  if (!self) {
    sn=which(network[,1]==network[,2])
    if (length(sn)>0)
      m=m[-sn,]
  }

  if (preorder) {
    # Ταξινόμηση κατά μειούμενο πλήθος στόχων - Προαιρετικό
    m1=matrix(ncol=2,nrow=0)
    colnames(m1)=colnames(m)

    # s: Ρυθμίζοντα γονίδια (source)
    s=names(sort(table(m[,1]),decreasing = T))

    for (i in seq_along(s)) {
      m1=rbind(m1,m[which(m[,1]==s[i]),])
    }
    m=m1
    rm(m1)
    # Τέλος ταξινόμησης
  } else {
    # s: Ρυθμίζοντα γονίδια (source)
    s=unique(m[,1])
  }

  # Συνολικά γονίδια που συμμετέχουν στα υποδίκτυα (είτε ρυθμίζοντα, είτε ρυθμιζόμενα)
  g=unique(c(m)) # ~ union(unique(m[,1]),unique(m[,2]))
  len_g=length(g)
  if (self)
    divl=len_g
  else
    divl=len_g-1

  # Γονίδια που είναι μόνο ρυθμιζόμενα (στόχοι dangling)
  d=setdiff(unique(m[,2]),s)

  if (verbose)
    print(sprintf("%d Genes - %d dangling",len_g,length(d)))

  # Πίνακας πιθανοτήτων μετάβασης (ρύθμισης εν προκειμένω)
  # H πηγή είναι η στήλη και οι στόχοι οι γραμμές,
  # συνεπώς κάθε στήλη έχει άθροισμα 1
  H=matrix(0,len_g,len_g)
  rownames(H)=colnames(H)=g

  # Τα dangling μοιράζονται ισοδύναμα σε όλα τα υπόλοιπα γονίδια
  dperc=1.0/divl
  H[g,d]=dperc
  if (!self) diag(H[d,d])=0

  # Ρυθμίζοντα γονίδια (source)
  # Αν θέλουμε μπορούμε να αποδώσουμε μια ελάχιστη αξία
  # στις υπόλοιπες, αδύνατες ή άγνωστες μεταβάσεις (ρυθμίσεις).
  # Πρέπει να είναι πολύ μικρότερη του ποσοστού των dangling (π.χ. 1/100)
  # παράμετροι letZeros, divider
  used_perc=0
  for (cur in s) {
    cur_targets=m[which(m[,1]==cur),2]
    if (!letZeros) {
        non_targets=setdiff(g,cur_targets)
        if (!self)
          non_targets=non_targets[non_targets!=cur] # Αφαιρώ το ίδιο το γονίδιο
        H[non_targets,cur]=dperc/divider
        used_perc=(length(non_targets)*(dperc/divider))
    }

    # Το υπόλοιπο μοιράζεται στις γνωστές μεταβάσεις (ρυθμίσεις)
    perc=(1.0-used_perc)/length(cur_targets)
    H[cur_targets,cur]=perc
  }

  # Εφαρμογή συντελεστή απόσβεσης
  H=damping*H+((1.0-damping)/divl)
  if (!self)
    diag(H)=0

  # Αρχικοποίηση [1,0,0,....]
  I=matrix(data=0,nrow=len_g,ncol=1)
  rownames(I)=rownames(H)
  colnames(I)=c("Rank")
  I[1,1]=1

  # Τεστ με πίνακα 8X8 για τον οποίο γνωρίζω τα τελικά αποτελέσματα
  # που είναι Ι=c(0.06,0,0675,0.03,0.0675,0.0975,0.2025,0.18,0.295)
  # len_g=8
  # H[1:8,1:8]=c(0,0.5,0.5,0,0,0,0,0,
  #     0,0,0,1,0,0,0,0,
  #     0,0.5,0,0,0.5,0,0,0,
  #     0,1/3,0,0,1/3,1/3,0,0,
  #     0,0,0,0,0,1/3,1/3,1/3,
  #     0,0,0,0,0,0,0,1,
  #     1/3,0,0,0,1/3,0,0,1/3,
  #     0,0,0,0,0,0.5,0.5,0)
  # H=H[1:8,1:8]
  # I[1:8]=1/8
  # I=I[1:8]

  # Power method
  for (i in seq_len(max_iterations)) {
    # if (i==60) {
    #   print("Break")
    # }
    I1=H%*%I
    e=sqrt(sum((I1 - I)^2))  # Ευκλίδεια απόσταση τρέχοντος βήματος από το προηγούμενο
    I=I1
    if (e<=threshold) break  # Ικανοποιητική σύγκλιση
  }
  if (verbose)
    print(sprintf("Finished in %d iterations",i))
  if (sorted)
    I=I[order(I,decreasing = T),,drop=F]
  return(I)
}
