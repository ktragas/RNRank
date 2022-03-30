#' @encoding UTF-8
#' Rank by gene importance
#'
#' @description Υπολογισμός σημαντικότητας γονιδίων ρυθμιστικού δικτύου με χρήση αλγορίθμου
#' παρόμοιου με τον PageRank της Google.
#'
#' @param srcm Πίνακας ονομάτων γονιδίων 2 στηλών. Κάθε γραμμή περιέχει ένα ρυθμίζον
#' κι ένα ρυθμιζόμενο γονίδιο.
#' @param max_iterations Ο μέγιστος αριθμός επαναλήψεων, αν δεν επιτευχθεί σύγκλιση νωρίτερα (minimum 10).
#' @param threshold Ευκλίδεια απόσταση μικρότερη ή ίση από αυτήν την τιμή μεταξύ των βημάτων,
#' θα τερματίσει την επανάληψη.
#' @param damping Συντελεστής απόσβεσης (0-1) που καθορίζει ποσοστό μη τυχαιότητας.
#' @param self Καθορίζει αν επιτρέπονται οι εγγραφές αυτορρύθμισης.
#' @param letZeros Καθορίζει αν επιτρέπουμε τυχαίες μεταβάσεις από τα ρυθμίζοντα γονίδια
#' προς άλλα που δεν καθορίζονται από τον πίνακά μας.
#' @param divider Καθορίζει διαιρέτη του ποσοστού που υπολογίζεται για τις μεταβάσεις
#' από τα μη ρυθμίζοντα γονίδια, ώστε το αποτέλεσμα της διαίρεσης να χρησιμοποιηθεί
#' για τις τυχαίες/άγνωστες μεταβάσεις από τα ρυθμίζοντα (minimum 10).
#' @param sorted Καθορίζει αν τα αποτελέσματα θα επιστραφούν ταξινομημένα κατά ποσοστό.
#' @param verbose Εμφανίζει κάποια μηνύματα κατά τη διαδικασία.
#' @param throwOnError Αν είναι FALSE, σε περίπτωση σοβαρού λάθους επιστρέφει NULL.
#' Αλλιώς καλείται η stop().
#'
#' @return Ονοματισμένος πίνακας ποσοστών (0-1) 1 στήλης.
#' @export
RNRank = function(srcm, damping=1.0, max_iterations=100, threshold=0,
                  self=F, letZeros=F, divider=100.0, sorted=T,
                  verbose=F, throwOnError=T)
{
  # Έλεγχος παραμέτρων
  if (!is.matrix(srcm) || ncol(srcm)<2 || nrow(srcm)<2 || typeof(srcm)!="character") {
    if (throwOnError)
      stop("Need character matrix with at least 2 rows and 2 columns. Can't continue")
    return(NULL)
  }

  warn=function(x) {
    msg=sprintf("Invalid value for \"%s\". Adjusted to %d",deparse(substitute(x)),x)
    warning(msg)
  }

  if (damping<0 || damping>1) {
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
  if (!self) {
    sn=which(srcm[,1]==srcm[,2])
    m=srcm[-sn,]
  } else
    m=srcm

  # Ταξινόμηση κατά μειούμενο πλήθος στόχων - Προαιρετικό
  m1=matrix(ncol=2,nrow=0)
  colnames(m1)=colnames(m)

  # s: Ρυθμίζοντα γονίδια (source)
  s=names(sort(table(m[,1]),decreasing = T))

  for (i in 1:length(s)) {
    m1=rbind(m1,m[which(m[,1]==s[i]),])
  }
  m=m1
  rm(m1)
  # Τέλος ταξινόμησης

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

  H=damping*H+((1.0-damping)/divl)
  if (!self)
    diag(H)=0

  # Αρχικοποίηση [1,0,0,....]
  I=matrix(data=0,nrow=len_g,ncol=1)
  rownames(I)=rownames(H)
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
  for (i in 1:max_iterations) {
    # if (i==60) {
    #   print("Break")
    # }
    I1=H%*%I
    e=sqrt(sum((I1 - I)^2))  # Ευκλίδεια απόσταση τρέχοντος βήματος από το προηγούμενο
    if (e<=threshold) break # Σύγκλιση
    I=I1
  }
  if (verbose)
    print(sprintf("Finished in %d iterations",i))
  if (sorted)
    I=I[order(I,decreasing = T),,drop=F]
  return(I)
}
