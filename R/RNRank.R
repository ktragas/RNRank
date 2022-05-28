#' RNRank: A package for gene ranking in inferred active regulatory network
#'
#' The RNRank package provides eight functions:
#' * [RNEAv3()],
#' * [RNRank()],
#' * [RNFfl()],
#' * [RNCircle()],
#' * [produce_reference()],
#' * [produce_reference_from_files()],
#' * [getGeneNames()] and
#' * [getPackageData()]
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
#' @details `favoredFFLs` and `favoredCircles` contain gene names and are treated
#' differently, depending on them being vectors or matrices and if any of them
#' is a matrix, depending on the number of columns. In any other case, they are ignored.
#'
#' * Vector or 1-column matrix: The regulations with any other gene are favored.
#' * 2 columns: Considered to contain directed edges. Exactly these regulations are favored.
#' * More columns: `favoredFFLs` is considered to contain FFLs in the first three columns
#'  (as the matrix returned from [RNFfl()]). The regulations C1->C2, C1->C3 and
#'  C2->C3 are favored. `favoredCircles` is considered to contain circles in the
#'  first three columns (as the matrix returned from [RNCircle()]).
#'  The regulations C1->C2, C2->C3 and C3->C1 are favored.
#'
#' @eval net_param()
#' @param max_iterations Maximum number of iterations, if not converging earlier (minimum 10).
#' @param threshold Euclidean distance between iterations less than or equal to this value
#' will terminate calculations
#' @param damping Damping factor (0-1) defining probability of non-randomness.
#' @param self Self regulations permitted (`TRUE`) or not (`FALSE`).
#' @param letZeros Set to `FALSE` to allow random regulations from regulating genes
#' to those not normally regulated.
#' @param divider Divides the probability used for dangling genes,
#' to use it for random regulations from regulating genes (minimum 10).
#' @param sorted Results sorted by descending probability (T) or unsorted (F).
#' @param preorder Sort input matrix, so that regulating genes with more targets are first.
#' @param favoredFFLs Regulations or genes, favored for better ranking. See 'Details'.
#' @param favoredCircles Regulations or genes, favored for better ranking. See 'Details'.
#' @param favorings Multiplier of favored regulations. Single value or numeric vector.
#' In case of vector, first value is for `favoredFFLs` and second for `favoredCircles`.
#' @param returned Function can return values, ranks or both.
#' @param ties.method	a character string specifying how ties are treated. See [base::rank()].
#' @eval common_params()
#'
#' @return If `returned` is "`values`", a named 1-column matrix of probabilities (0-1) is returned.
#' If `returned` is "`ranks`", a named 1-column matrix of ranks (highest value->first ranking) is returned.
#' If `returned` is "`both`", a list containing the two matrices is returned.
#' @export
#'
#' @examples
RNRank = function(network, damping=0.85, max_iterations=100, threshold=0,
                  self=T, letZeros=F, divider=100.0, sorted=T, preorder=F,
                  favoredFFLs=NULL, favoredCircles=NULL, favorings=1.1,
                  returned=c("values","ranks","both"),
                  ties.method = "min", verbose=F, throwOnError=T)
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
      m1=rbind(m1,m[m[,1]==s[i],])
    }
    m=m1
    rm(m1)
    # Τέλος ταξινόμησης
  } else {
    # s: Ρυθμίζοντα γονίδια (source)
    s=unique(m[,1])
  }

  # Συνολικά γονίδια που συμμετέχουν στο υποδίκτυο (είτε ρυθμίζοντα, είτε ρυθμιζόμενα)
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
    cur_targets=m[m[,1]==cur,2]
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

  # Πριμοδοτήσεις (των FFL και των κύκλων π.χ.)
  if (is.numeric(favorings)) {
    if (length(favorings)==1)
       favorings=rep(favorings,2)
    for (i in 1:2) {
       if (favorings[i]<=0 || favorings[i]==1) next
       if (i==1) {
          f=favoredFFLs
          what="ffl"
       } else {
          f=favoredCircles
          what="circle"
       }
       if (!is.null(f))
          H=adjustTransitions(f,H,favorings[i],what)
    }
  }

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

  returned=match.arg(returned)
  if (returned!="values")
    R=rank(-I,ties.method = ties.method)

  if (returned=="ranks") {
    return(R)
  } else if (returned=="both") {
    return(list(values=I,ranks=R))
  } else # "values"
    return(I)
}
