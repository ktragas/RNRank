#' @title
#' Regulatory Network Enrichment Analysis (v3)
#'
#' @description
#' Infers active regulatory networks from gene expression data using a combination
#' of prior knowledge and enrichment analysis
#'
#' @details
#' Transcription factors (TF) control a multitude of constitutive, cell-specific,
#' developmental, proliferative or homeostatic processes in the cells of all known
#' organisms. Due to their central role in gene regulation a considerable number
#' of human diseases have been associated with TF function, including metabolic,
#' autoimmune disorders and cancer. Despite the progress over the last years,
#' towards the identification of TF gene targets, either with experimental or
#' in silico approaches, we are still far from accurately reconstructing the
#' hierarchy of transcriptional regulators from genome-wide data. In this work
#' we are trying to overcome the existing limitations for the inference and study
#' of regulatory networks of TF-interactions and their subsequent use in
#' enriching our understanding of key biological processes. Given a differential
#' expression file RNEA extracts lists of prioritized TFs, and possibly miRNAs.
#' Most importantly a regulatory subnetwork is extracted showing how the activated
#' regulators interact with their target genes and between each other.
#' In this way, given a genome-wide expression experiment, RNEA prioritizes
#' important regulatory components.
#' It also succeeds to reconstruct meaningful subnetworks of gene regulation,
#' offering further ways of analyzing the data based on network theory.
#' A detailed guide describing RNEA's use can be found [here](https://sites.google.com/a/fleming.gr/rnea/manual).
#'
#' @param filename Differential expression file (Cuffdiff, 3-column, 1-column)
#' @param species Organism. "Human","Mouse" datasets exist in package.
#' If any other, the relevant reference file must be provided in `reference_dir`.
#' @param internal_data If `TRUE` the package's internal datasets are used
#' @param miRNA_data Use miRNA data
#' @param FC_threshold log2FC threshold
#' @param PV_threshold p-value threshold
#' @param output_dir Directory to save output to
#' @param output Prefix for output files
#' @param reference_dir Directory containing reference files for unsupported species
#' @param ffl Locate and save Feed-Forward Loops by calling [RNFfl()]
#' @param circles Locate and save circles by calling [RNCircle()]
#' @param rank Call RNRank() on the resulting network and save an extra file with ranks
#' @param favorFFL If ranking, favor genes participating in FFLs
#' @param favorCircles If ranking, favor genes participating in circles
#' @param ... Arguments to be passed to [RNRank()]
#' @eval common_params()
#'
#' @return "Analysis completed!"
#' @export
#'
#' @examples
RNEAv3<-function(filename,species,internal_data=T, miRNA_data=F,
                 FC_threshold=1,PV_threshold=0.05,
                 output_dir=".",output="Output",
                 reference_dir="ReferenceFiles",ffl=T,circles=T,rank=T,
                 favorFFL=T, favorCircles=T, ..., verbose=F, throwOnError=T){

  # Κώστας - Ο έλεγχος των παραμέτρων γίνεται άμεσα
  if (internal_data==TRUE && (species %in% c("Mouse","Human"))==FALSE)
    warning("\"Mouse\" and \"Human\" are the only species supported internally as of now");

  # Μετά από τον έλεγχο παραμέτρων, όλα τα warnings αγνοούνται μέχρι το τέλος της συνάρτησης
  withr::local_options(list(warn = -1))

  # Αν δεν υπάρχει το directory εξόδου, το δημιουργώ
  if (dir.exists(output_dir)==FALSE)
    dir.create(output_dir)
  # Αν προϋπήρχε ή πέτυχε η δημιουργία το προσθέτω στο πρόθεμα
  if (dir.exists(output_dir)==TRUE)
    output=file.path(output_dir,output)

  if (internal_data && species %in% c("Human","Mouse")) {
    # Αν έχουν ζητηθεί τα δεδομένα που περιέχονται στο πακέτο για άνθρωπο ή ποντίκι
    if (species=="Human")
      TF=TF_human
     else   # species=="Mouse"
      TF=TF_mouse

    if (miRNA_data) {
      if (species=="Human")
        miRNA=miRNA_human
       else   # species=="Mouse"
        miRNA=miRNA_mouse
    }
  } else { # Don't use internal human or mouse data
   	# TF_ref=file.path(reference_dir,paste(species,"_TF_plus_oreganno_trrust.tsv",sep=""));
   	TF_ref=file.path(reference_dir,paste0(species,"_TF_Reference.tsv"));
   	TF=read.table(TF_ref,sep="\t",fill=T);
  	rownames(TF)=TF[,2];
  	if (miRNA_data) {
    	miRNA_ref=file.path(reference_dir,paste0(species,"_miRNAs.tsv"));
    	miRNA=read.table(miRNA_ref,sep="\t",fill=T);
    	rownames(miRNA)=miRNA[,2];
  	}
  }
  TF_genes=TF[,2]

	Input=read.table(filename,sep="\t",header=T);
	ncolInput=ncol(Input);
	if(ncolInput==14) {
		Total_number_of_genes=length(unique(Input[,3]));
	} else {
	  # Για όλες τις άλλες περιπτώσεις θεωρούμε ότι τα γονίδια είναι στην πρώτη στήλη
	  Total_number_of_genes=length(unique(Input[,1]));
	}

	TF_counts=cbind(TF[,1],Total_number_of_genes,0,0,0,0);
	colnames(TF_counts)=c("Targets","Genes","IsDE?","TargetsDE","TargetsUp","TargetsDown");
	rownames(TF_counts)=TF[,2];
	if (miRNA_data) {
  	miRNA_counts=cbind(miRNA[,1],Total_number_of_genes,0,0,0);
  	colnames(miRNA_counts)=c("Targets","Genes","TargetsDE","TargetsUp","TargetsDown");
  	rownames(miRNA_counts)=miRNA[,2];
	}
	Number_of_DE_genes=0;
	Number_of_Up_genes=0;
	Number_of_Down_genes=0;
	DE_genes=matrix(data=0,nrow=1);
	UP_genes=matrix(data=0,nrow=1);
	DOWN_genes=matrix(data=0,nrow=1);
	# Network=matrix(ncol=2,data=""); # Κώστας - Δημιουργούσε μια κενή γραμμή στα αποτελέσματα
	Network=matrix(ncol=2,nrow=0);
	colnames(Network)=c("Source","Target");

	# Κώστας (start)
	if (ncolInput==14) { g_idx=3; fc_idx=10; pv_idx=12; }
	else {
	  g_idx=1;
	  if (ncolInput==3) { fc_idx=2; pv_idx=3; }
	}

	for (line in seq_len(nrow(Input))){ # Στον παλιό κώδικα έλεγε for(line in 2:nrow(Input)) - πρώην header???
	  PV=Input[line,pv_idx]
	  FC=Input[line,fc_idx]
	  # sanity checks
	  if (!is.numeric(PV) || !is.finite(PV)) next;
	  if (!is.numeric(FC) || !is.finite(FC)) next;
	  # threshold checks
	  if (PV>PV_threshold) next;
	  if (FC>-FC_threshold && FC<FC_threshold) next;

	  # Εδώ έρχονται μόνο DE
	  gene=as.character(Input[line, g_idx])
	  if (length(gene)<=0)
	    gene="?";

	  Number_of_DE_genes=Number_of_DE_genes+1;
	  DE_genes[Number_of_DE_genes]=gene;

	  # Αν το input έχει μια στήλη, θεωρούνται upregulated by default
	  if (ncolInput==1 || Input[line,fc_idx]>=FC_threshold) {
	    Number_of_Up_genes=Number_of_Up_genes+1;
	    UP_genes[Number_of_Up_genes]=gene;
	  } else { # ncolInput!=1 && Input[line,fc_idx]<=-FC_threshold
	    Number_of_Down_genes=Number_of_Down_genes+1;
	    DOWN_genes[Number_of_Down_genes]=gene;
	  }
	}
	# Κώστας (end)

	DE_genes=unique(DE_genes);
	UP_genes=unique(UP_genes);
	DOWN_genes=unique(DOWN_genes);
	Number_of_DE_genes=length(DE_genes);
	Number_of_Up_genes=length(UP_genes);
	Number_of_Down_genes=length(DOWN_genes);
	if (verbose) {
  	print(paste(Number_of_DE_genes,"differentially expressed genes"));
  	print(paste(Number_of_Up_genes,"upregulated"));
  	print(paste(Number_of_Down_genes,"downregulated"));
	}

	# Kώστας (start)
	Genes=list(up=UP_genes,down=DOWN_genes)
	IsDE=list(up=1,down=-1)
	if (verbose) count=0
	for (direction in c("up","down")) {
	  # zar : Διαφορικά εκφρασμένα ρυθμιστικά γονίδια
	  # (δηλ. που βρίσκονται στη στήλη 2 των αρχείων reference)
	  zar = intersect(TF_genes, Genes[[direction]])
    for (z in zar) { # z : parent (=1ου επιπέδου) γονίδιο
      if (verbose) {
        count=count+1
        cat(sprintf("[%d] %s\r",count,stringr::str_pad(z,15,"right")))
      }

      # Κάθε εγγραφή περιέχει διαφορετικό πλήθος ρυθμιζόμενων γονιδίων,
      # που καθορίζεται στη στήλη 1 των αρχείων reference
      # και βρίσκονται από την στήλη 3 και πέρα
      last_index = 2 + TF[z, 1]
      regulated_genes = TF[z, 3:last_index]
      deup = intersect(regulated_genes, UP_genes)
      dedown = intersect(regulated_genes, DOWN_genes)
      deup_len = length(deup)
      dedown_len = length(dedown)
      # Parent - children
      if (deup_len > 0 || dedown_len > 0) {
        Network = rbind(Network,
                  if_else(deup_len>0,cbind(z,deup),NULL),
                  if_else(dedown_len>0,cbind(z,dedown),NULL))
      }

      TF_counts[z, 3] = IsDE[[direction]]
      TF_counts[z, 4] = TF_counts[z, 4] + deup_len + dedown_len
      TF_counts[z, 5] = TF_counts[z, 5] + deup_len
      TF_counts[z, 6] = TF_counts[z, 6] + dedown_len

      child_TFs=intersect(TF[z,3:last_index],TF_genes)
      # child_TF_indices=match(child_TFs,TF[z,3:last_index])+2

      # Για κάθε ρυθμιζόμενο γονίδιο που ρυθμίζει κι αυτό κάποια γονίδια (είναι δηλαδή TF)
      for (gene in child_TFs) {
        inner_last_index=2 + TF[gene, 1]
        deup = intersect(TF[gene, 3:inner_last_index], UP_genes)
        dedown = intersect(TF[gene, 3:inner_last_index], DOWN_genes)
        deup_len = length(deup)
        dedown_len = length(dedown)
        # Parent - Child & Child - GrandChildren
        # Αν το ρυθμιζόμενο γονίδιο είναι DE, ίσως έχει ήδη προστεθεί
        # στο προηγούμενο στάδιο, αλλά είναι πιο εύκολο να αφαιρεθεί
        # η διπλή εγγραφή με unique() στο τέλος, από το να γίνονται
        # έλεγχοι εδώ
        if (deup_len > 0 || dedown_len > 0) {
          Network = rbind(Network, c(z,gene),
            if_else(deup_len>0,cbind(gene,deup),NULL),
            if_else(dedown_len>0,cbind(gene,dedown),NULL))
        }

        TF_counts[z, 1] = TF_counts[z, 1] + TF[gene, 1]
        TF_counts[z, 4] = TF_counts[z, 4] + deup_len + dedown_len
        TF_counts[z, 5] = TF_counts[z, 5] + deup_len
        TF_counts[z, 6] = TF_counts[z, 6] + dedown_len
      } # for (gene...
    } # for (z...
	} # for (direction...

	if (miRNA_data) {
  	for(mir in seq_len(nrow(miRNA))){
  	  z=rownames(miRNA)[mir];
  	  last_index=2+miRNA[z,1];
  	  regulated_genes=miRNA[z,3:last_index]
  	  deup=intersect(regulated_genes,UP_genes);
  	  dedown=intersect(regulated_genes,DOWN_genes);
  	  de=intersect(regulated_genes,DE_genes);
  	  if(length(de)>0)
  	        Network=rbind(Network,cbind(z,de));
      # Πάει μόνο σε ένα επίπεδο, γιατί αν κάποιο από τα διαφορικά εκφρασμένα
  	  # ρυθμιζόμενα γονίδια ήταν ρυθμιστής, θα είχαμε ήδη ασχοληθεί μαζί του
  	  # στο προηγούμενο στάδιο
  	  miRNA_counts[z,3]=miRNA_counts[z,3]+length(de);
  	  miRNA_counts[z,4]=miRNA_counts[z,4]+length(deup);
  	  miRNA_counts[z,5]=miRNA_counts[z,5]+length(dedown);
  	} # for (mir...
	}

	DE_rows=which(TF_counts[,3]!=0); # Διαφορικά εκφρασμένα γονίδια
	if(length(DE_rows)==0){
	  print("No TF's targets found deregulated")
	} else {
	  # print(TF_counts[which(TF_counts[,3]!=0),]);
	  ResultsTF=matrix(nrow=nrow(TF_counts[DE_rows,]),ncol=3);
	  rownames(ResultsTF)=rownames(TF_counts)[DE_rows];
	  colnames(ResultsTF)=c("DE_qvalue","UP_qvalue","DOWN_qvalue");
	  ResultsTFtemp=matrix(nrow=nrow(TF_counts[DE_rows,]),ncol=3);
	  rownames(ResultsTFtemp)=rownames(TF_counts)[DE_rows];
	  colnames(ResultsTFtemp)=c("DE_pvalue","UP_pvalue","DOWN_pvalue");
	  for(i in 1:length(DE_rows)){
	    r=DE_rows[i];
	    # Hypergeometric Distr.          # White_Drawn   Total_White                 Total_Black                  Total_Drawn
	                                     # TargetsDE     TotalDE                     TotalNDE                     Targets
	    ResultsTFtemp[i,1]<-stats::phyper(TF_counts[r,4],Number_of_DE_genes,(TF_counts[r,2]-Number_of_DE_genes), TF_counts[r,1],lower.tail=FALSE);
	    ResultsTFtemp[i,2]<-stats::phyper(TF_counts[r,5],Number_of_Up_genes,(TF_counts[r,2]-Number_of_Up_genes), TF_counts[r,1],lower.tail=FALSE);
	    ResultsTFtemp[i,3]<-stats::phyper(TF_counts[r,6],Number_of_Down_genes,(TF_counts[r,2]-Number_of_Down_genes), TF_counts[r,1],lower.tail=FALSE);
	  }
	  ResultsTF[,1]=stats::p.adjust(ResultsTFtemp[,1],method="fdr");
	  ResultsTF[,2]=stats::p.adjust(ResultsTFtemp[,2],method="fdr");
	  ResultsTF[,3]=stats::p.adjust(ResultsTFtemp[,3],method="fdr");
	}

	if (miRNA_data) {
  	if(length(which(miRNA_counts[,3]!=0))==0){
  		print("No miRNA's targets found deregulated")
  	} else {
  # 		print(miRNA_counts[which(miRNA_counts[,3]!=0),]);
  		ResultsmiRNA=matrix(nrow=nrow(miRNA_counts[which(miRNA_counts[,3]!=0),]),ncol=3);
  		rownames(ResultsmiRNA)=rownames(miRNA_counts)[which(miRNA_counts[,3]!=0)];
  		colnames(ResultsmiRNA)=c("DE_qvalue","UP_qvalue","DOWN_qvalue");
  		ResultsmiRNAtemp=matrix(nrow=nrow(miRNA_counts[which(miRNA_counts[,3]!=0),]),ncol=3);
  		rownames(ResultsmiRNAtemp)=rownames(miRNA_counts)[which(miRNA_counts[,3]!=0)];
  		colnames(ResultsmiRNAtemp)=c("DE_pvalue","UP_pvalue","DOWN_pvalue");
  		for(i in 1:length(which(miRNA_counts[,3]!=0))){
  			r=which(miRNA_counts[,3]!=0)[i];
  			ResultsmiRNAtemp[i,1]<-stats::phyper(miRNA_counts[r,3],Number_of_DE_genes,(miRNA_counts[r,2]-Number_of_DE_genes), miRNA_counts[r,1],lower.tail=FALSE);
  			ResultsmiRNAtemp[i,2]<-stats::phyper(miRNA_counts[r,4],Number_of_Up_genes,(miRNA_counts[r,2]-Number_of_Up_genes), miRNA_counts[r,1],lower.tail=FALSE);
  			ResultsmiRNAtemp[i,3]<-stats::phyper(miRNA_counts[r,5],Number_of_Down_genes,(miRNA_counts[r,2]-Number_of_Down_genes), miRNA_counts[r,1],lower.tail=FALSE);
  		}
  		ResultsmiRNA[,1]=stats::p.adjust(ResultsmiRNAtemp[,1],method="fdr");
  		ResultsmiRNA[,2]=stats::p.adjust(ResultsmiRNAtemp[,2],method="fdr");
  		ResultsmiRNA[,3]=stats::p.adjust(ResultsmiRNAtemp[,3],method="fdr");
  	}
	}
	Network=unique(Network)

	write.csv(Network,file=paste0(output,"_Network.csv"),quote=F,row.names=F);
	if (ffl) {
	  outffl=RNFfl(Network, verbose, throwOnError)
	  if (nrow(outffl)>0)
	    write.csv(outffl,file=paste0(output,"_FFLs.csv"),quote = F,row.names = F)
	}
	if (circles) {
	  outcircle=RNCircle(Network,verbose,throwOnError)
	  if (nrow(outffl)>0)
	    write.csv(outcircle,file=paste0(output,"_Circles.csv"),quote = F,row.names = F)
  }

	if (rank) {
	  P=RNRank(Network,...,favoredFFLs = ifelse(favorFFL & ffl,outffl,NULL),
	           favoredCircles=ifelse(favorCircles & circles,outcircle,NULL),
	           verbose=verbose,throwOnError = throwOnError)
	  write.csv(P,file=paste0(output,"_Ranks.csv"),quote = F,row.names = T)
	}

	tempnames=c("Regulator",colnames(ResultsTF));
	ResultsTF=cbind(rownames(ResultsTF),as.data.frame(ResultsTF));
	colnames(ResultsTF)=tempnames;
	write.csv(as.data.frame(ResultsTF),file=paste0(output,"_TF_Enrichment.csv"),row.names=F);

	if (miRNA_data) {
  	tempnames=c("miRNA",colnames(ResultsmiRNA));
  	ResultsmiRNA=cbind(rownames(ResultsmiRNA),as.data.frame(ResultsmiRNA));
  	colnames(ResultsmiRNA)=tempnames;
  	write.csv(as.data.frame(ResultsmiRNA),file=paste0(output,"_miRNA_Enrichment.csv"),row.names=F);
  }

	# Return
	print("Analysis completed!");
}



