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



# Προσαρμογή μεταβάσεων σύμφωνα με το δίκτυο
#
# m: Ακμές ή γονίδια που θα πριμοδοτηθούν.
#  Πίνακας με μία στήλη ή διάνυσμα: θα οδηγήσει σε πριμοδότηση κάθε σύνδεσης
#  των στοιχείων (Η[C1,] και H[,C1])
#  Πίνακας με δύο στήλες: θεωρώ ότι περιέχει κατευθυνόμενες ακμές
#  Πριμοδοτούνται ακριβώς αυτές οι μεταβάσεις (Η[C2,C1])
#  Πίνακας με τρεις ή παραπάνω στήλες: Κάθε γραμμή του τη μετατρέπω σε 3 γραμμές
#  2 στηλών, δηλαδή στις αντίστοιχες ακμές. Aν η παράμετρος func ορίζει ότι στις 3
#  πρώτες στήλες περιέχει FFL, τελικά θα πριμοδοτηθούν οι μεταβάσεις H[C2,C1],
#  H[C3,C2] και H[C3,C1]. Αν περιέχει κύκλους, τότε θα πριμοδοτηθούν οι μεταβάσεις
#  H[C2,C1], H[C3,C2] και H[C1,C3].
# Η: Πίνακας πιθανοτήτων μετάβασης όπως έχει ετοιμαστεί από την RNRank().
# w: Συντελεστής με τον οποίο πολλαπλασιάζονται τα κελιά που πρέπει
#  >1 Πριμοδότηση, <1 Ποινή
# func: Καθορίζει αν ο πίνακας m τριών στηλών και άνω αναπαριστά FFL ή κύκλο
adjustTransitions<-function(m,H,w=1,func=c("ffl","circle"))
{
   if (is.vector(m))
      m=matrix(unique(m),ncol=1)
   if (!is.matrix(m))
      return(H)
   nc=ncol(m)
   if (nc<1)
      return(H)
   if (nc>=3) {
      func=match.arg(func)
      if (func=="ffl")
         m=unique(rbind(m[,1:2],m[,c(1,3)],m[,2:3]))
      else #circle
         m=unique(rbind(m[,1:2],m[,2:3],m[,c(3,1)]))
      nc=2
   }

   if (nc==2) {
      for (r in seq_len(nrow(m))) {
         s=m[r,1]
         t=m[r,2]
         H[t,s]=H[t,s]*w
      }
   } else if (nc==1) {
      i=which(rownames(H) %in% m[,1])
      h=H[-i,-i] # Κρατάω τα στοιχεία που δεν πρέπει να αλλάξουν
      H=H*w      # Πολλαπλασιάζω όλον τον πίνακα
      H[-i,-i]=h # Επαναφέρω αυτά που κράτησα
      # Αν το έκανα με κάποιον άλλον vectorized τρόπο κάποια στοιχεία
      # θα πολλαπλασιάζονταν πάνω από μια φορά, ή κάποια δεν θα
      # πολλαπλασιάζονταν ενώ έπρεπε
   }

   # Κλιμάκωση των τιμών ώστε κάθε στήλη να δίνει άθροισμα 1
   H=apply(H,2,function(x) x/sum(x))
   return(H)
}

# Hack για να επιστρέφεται και NULL χωρίς λάθος
if_else<-function(s,a,b) switch((s)+1,b,a)

