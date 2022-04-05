#' @title
#' Regulatory Network Enrichment Analysis (v3)
#'
#' @description
#' Infers active regulatory networks from gene expression data using a combination
#' of prior knowledge and enrichment analysis
#'
#' @details
#' Transcription factors (TF) control a multitude of constitutive, cell-specific, developmental, proliferative or homeostatic processes in the cells of all known organisms. Due to their central role in gene regulation a considerable number of human diseases have been associated with TF function, including metabolic, autoimmune disorders and cancer. Despite the progress over the last years, towards the identification of TF gene targets, either with experimental or in silico approaches, we are still far from accurately reconstructing the hierarchy of transcriptional regulators from genome-wide data. In this work we are trying to overcome the existing limitations for the inference and study of regulatory networks of TF-interactions and their subsequent use in enriching our understanding of key biological processes. Furthermore we are combining our approach with state of the art functional enrichment analyses in order to create a tool, called Regulatory Network Enrichment Analysis (RNEA) that will prioritize transcriptional and functional characteristics of a genome-wide expression experiment.
#'
#' Given a differential expression file RNEA extracts lists of prioritized TFs, miRNAs, KEGG pathways and categories and GO terms.
#' Most importantly a regulatory subnetwork is extracted showing how the activated regulators interact with their target genes and between each other.
#' In this way, given a genome-wide expression experiment, RNEA prioritizes important regulatory and functional components.
#' It also succeeds to reconstruct meaningful subnetworks of gene regulation, offering further ways of analyzing the data based on network theory.
#' A detailed guide describing RNEA's use can be found [here](https://sites.google.com/a/fleming.gr/rnea/manual).
#'  Additionally RNEA can extract a global regularory/functional network with the "activated" regulators, their target genes and functional categories whose members are overrepresented among differentially expressed genes.
#'  This network, although sometimes huge, gives a complete view in both functional and regulatory components affected at the system studied.
#'
#' @param filename Differential expression file (Cuffdiff, 3-column, 1-column)
#' @param identifier default "GeneName" or (not working) "RefSeq"
#' @param species Organism: "Human" or "Mouse"
#' @param FC_threshold log2FC threshold
#' @param PV_threshold p-value threshold
#' @param network Needed output: "global" or (default) "regulatory" or "functional"
#' @param output_dir Directory to save output to
#' @param output Prefix for output files
#' @param type_of_output Type of output: (default) "csv" or "html"
#' @param reference_dir Directory containing reference files for unsupported species
#' @param rank Call RNRank() on the resulting network and save an extra file with ranks
#' @param ... Arguments to be passed to [RNRank()]
#' @return "Analysis completed!"
#' @export
#'
#' @examples
RNEAv3<-function(filename,identifier="GeneName",species,
                 FC_threshold=1,PV_threshold=0.05,network="regulatory",
                 output_dir=".",output="Output",type_of_output="csv",
                 reference_dir="ReferenceFiles",rank=T, ...){

  # Κώστας - Ο έλεγχος των παραμέτρων να γίνεται άμεσα
  if ((identifier %in% c("GeneName","Refseq"))==FALSE)
    warning("Please select a supported gene identifier (\"GeneName\" or \"Refseq\")");
  if ((species %in% c("Mouse","Human"))==FALSE)
    warning("\"Mouse\" and \"Human\" are the only species supported as of now");
  if ((type_of_output %in% c("html","csv"))==FALSE)
    warning("Please select type of output to be either \"html\" or \"csv\"")
  if ((network %in% c("global","functional","regulatory"))==FALSE)
    warning("Please select network to be \"global\",\"functional\" or \"regulatory\"")
  needFunctional=network %in% c("global","functional")
  needRegulatory=network %in% c("global","regulatory")

  # Μετά από τον έλεγχο παραμέτρων, όλα τα warnings αγνοούνται
  options(warn=-1);

  if(identifier=="Refseq"){
		Refseq=utils::read.table("ReferenceFiles/Refseq2Gene.txt",row.names=2,sep="\t");
  }

  if (species %in% c("Human","Mouse")) {
    # Για άνθρωπο και ποντίκι τα δεδομένα περιέχονται στο πακέτο
    if (species=="Human") {
#      utils::data(TF_human,miRNA_human,kegg_human,keggcat_human,gos_human)
      TF=TF_human
      miRNA=miRNA_human
      kegg=kegg_human
      keggcat=keggcat_human
      gos=gos_human
    } else {  # species=="Mouse"
#      utils::data(TF_mouse,miRNA_mouse,kegg_mouse,keggcat_mouse,gos_mouse)
      TF=TF_mouse
      miRNA=miRNA_mouse
      kegg=kegg_mouse
      keggcat=keggcat_mouse
      gos=gos_mouse
    }
  } else { # Not human or mouse
   	TF_ref=file.path(reference_dir,paste(species,"_TF_plus_oreganno_trrust.tsv",sep=""));
 	  miRNA_ref=file.path(reference_dir,paste(species,"_miRNAs.tsv",sep=""));
 	  kegg_ref=file.path(reference_dir,paste(species,"_kegg.tsv",sep=""));
 	  keggcat_ref=file.path(reference_dir,paste(species,"_keggcat.tsv",sep=""));
   	GOs_ref=file.path(reference_dir,paste(species,"_GOs.tsv",sep=""));

   	TF=utils::read.table(TF_ref,sep="\t",fill=T);
  	rownames(TF)=TF[,2];
  	miRNA=utils::read.table(miRNA_ref,sep="\t",fill=T);
  	rownames(miRNA)=miRNA[,2];
  	kegg=utils::read.table(kegg_ref,sep="\t",fill=T);
  	rownames(kegg)=kegg[,2];
  	keggcat=utils::read.table(keggcat_ref,sep="\t",fill=T);
  	rownames(keggcat)=keggcat[,2];
   	gos=utils::read.table(GOs_ref,sep="\t",fill=T);
   	rownames(gos)=gos[,2];
  }

	Input=utils::read.table(filename,sep="\t",header=T);
	if(ncol(Input)==1){
	  Total_number_of_genes=length(unique(Input[,1]));
	}
	if(ncol(Input)==3){
		Total_number_of_genes=length(unique(Input[,1]));
	}else if(ncol(Input)==14){
		Total_number_of_genes=length(unique(Input[,3]));
	}

	TF_counts=cbind(TF[,1],Total_number_of_genes,0,0,0,0);
	colnames(TF_counts)=c("Targets","Genes","IsDE?","TargetsDE","TargetsUp","TargetsDown");
	rownames(TF_counts)=TF[,2];
	miRNA_counts=cbind(miRNA[,1],Total_number_of_genes,0,0,0);
	colnames(miRNA_counts)=c("Targets","Genes","TargetsDE","TargetsUp","TargetsDown");
	rownames(miRNA_counts)=miRNA[,2];
	kegg_counts=cbind(kegg[,1],Total_number_of_genes,0,0,0);
	colnames(kegg_counts)=c("Members","Genes","MembersDE","MembersUp","MembersDown");
	rownames(kegg_counts)=kegg[,2];
	keggcat_counts=cbind(keggcat[,1],Total_number_of_genes,0,0,0);
	colnames(keggcat_counts)=c("Members","Genes","MembersDE","MembersUp","MembersDown");
	rownames(keggcat_counts)=keggcat[,2];
 	go_counts=cbind(gos[,1],Total_number_of_genes,0,0,0);
 	colnames(go_counts)=c("Members","Genes","MembersDE","MembersUp","MembersDown");
 	rownames(go_counts)=gos[,2];
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
	ncolInput=ncol(Input);
	if (ncolInput==14) { g_idx=3; fc_idx=10; pv_idx=12; }
	else if (ncolInput==3) { g_idx=1; fc_idx=2; pv_idx=3; }
	else if (ncolInput==1) { g_idx=1; }

	for (line in 1:nrow(Input)){ # Στον παλιό κώδικα έλεγε for(line in 2:nrow(Input)) - πρώην header???
	  if (Input[line,pv_idx]>PV_threshold) next;
	  if (Input[line,fc_idx]>-FC_threshold && Input[line,fc_idx]<FC_threshold) next;

	  # Εδώ έρχονται μόνο DE
	  if (identifier=="GeneName")
	    gene=as.character(Input[line, g_idx])
	  else if (identifier=="Refseq")
	    gene=as.character(Refseq[Input[line,g_idx],1])
	  else
	    gene="?";
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

	# if(ncol(Input)==14){
	#   for(line in 2:nrow(Input)){
	#     if((Input[line,12]<=PV_threshold)&&(Input[line,10]>=FC_threshold)){
	#       Number_of_DE_genes=Number_of_DE_genes+1;
	#       Number_of_Up_genes=Number_of_Up_genes+1;
	#       if(identifier=="GeneName"){
	#         DE_genes[Number_of_DE_genes]=as.character(Input[line,3]);
	#         UP_genes[Number_of_Up_genes]=as.character(Input[line,3]);
	#       }else if(identifier=="Refseq"){
	#         DE_genes[Number_of_DE_genes]=as.character(Refseq[Input[line,3],1]);
	#         UP_genes[Number_of_Up_genes]=as.character(Refseq[Input[line,3],1]);
	#       }else{
	#         print("Please select a supported gene identifier (\"GeneName\" or \"Refseq\")");
	#       }
	#     }
	#     if((Input[line,12]<=PV_threshold)&&(Input[line,10]<=-FC_threshold)){
	#       Number_of_DE_genes=Number_of_DE_genes+1;
	#       Number_of_Down_genes=Number_of_Down_genes+1;
	#       if(identifier=="GeneName"){
	#         DE_genes[Number_of_DE_genes]=as.character(Input[line,3]);
	#         DOWN_genes[Number_of_Down_genes]=as.character(Input[line,3]);
	#       }else if(identifier=="Refseq"){
	#         DE_genes[Number_of_DE_genes]=as.character(Refseq[Input[line,3],1]);
	#         DOWN_genes[Number_of_Down_genes]=as.character(Refseq[Input[line,3],1]);
	#       }
	#       else{
	#         print("Please select a supported gene identifier (\"GeneName\" or \"Refseq\")");
	#       }
	#     }
	#   }
	# }
	# if(ncol(Input)==3){
	#   for(line in 2:nrow(Input)){
	#     if((Input[line,3]<=PV_threshold)&&(Input[line,2]>=FC_threshold)){
	#       Number_of_DE_genes=Number_of_DE_genes+1;
	#       Number_of_Up_genes=Number_of_Up_genes+1;
	#       if(identifier=="GeneName"){
	#         DE_genes[Number_of_DE_genes]=as.character(Input[line,1]);
	#         UP_genes[Number_of_Up_genes]=as.character(Input[line,1]);
	#       }else if(identifier=="Refseq"){
	#         DE_genes[Number_of_DE_genes]=as.character(Refseq[Input[line,1],1]);
	#         UP_genes[Number_of_Up_genes]=as.character(Refseq[Input[line,1],1]);
	#       }else{
	#         print("Please select a supported gene identifier (\"GeneName\" or \"Refseq\")");
	#       }
	#     }
	#     if((Input[line,3]<=PV_threshold)&&(Input[line,2]<=-FC_threshold)){
	#       Number_of_DE_genes=Number_of_DE_genes+1;
	#       Number_of_Down_genes=Number_of_Down_genes+1;
	#       if(identifier=="GeneName"){
	#         DE_genes[Number_of_DE_genes]=as.character(Input[line,1]);
	#         DOWN_genes[Number_of_Down_genes]=as.character(Input[line,1]);
	#       }else if(identifier=="Refseq"){
	#         DE_genes[Number_of_DE_genes]=as.character(Refseq[Input[line,1],1]);
	#         DOWN_genes[Number_of_Down_genes]=as.character(Refseq[Input[line,1],1]);
	#       }
	#       else{
	#         print("Please select a supported gene identifier (\"GeneName\" or \"Refseq\")");
	#       }
	#     }
	#   }
	# }
	#
	# if(ncol(Input)==1){
	#   for(line in 2:nrow(Input)){
	#     Number_of_DE_genes=Number_of_DE_genes+1;
	#     Number_of_Up_genes=Number_of_Up_genes+1;
	#     if(identifier=="GeneName"){
	#       DE_genes[Number_of_DE_genes]=as.character(Input[line,1]);
	#       UP_genes[Number_of_Up_genes]=as.character(Input[line,1]);
	#     }else if(identifier=="Refseq"){
	#       DE_genes[Number_of_DE_genes]=as.character(Refseq[Input[line,1],1]);
	#       UP_genes[Number_of_Up_genes]=as.character(Refseq[Input[line,1],1]);
	#     }else{
	#       print("Please select a supported gene identifier (\"GeneName\" or \"Refseq\")");
	#     }
	#   }
	# }
	#

	# Κώστας (end)

	DE_genes=unique(DE_genes);
	UP_genes=unique(UP_genes);
	DOWN_genes=unique(DOWN_genes);
	Number_of_DE_genes=length(DE_genes);
	Number_of_Up_genes=length(UP_genes);
	Number_of_Down_genes=length(DOWN_genes);
	print(paste(Number_of_DE_genes,"differentially expressed genes"));
	print(paste(Number_of_Up_genes,"upregulated"));
	print(paste(Number_of_Down_genes,"downregulated"));

	# Kώστας (start)
	Genes=list(up=UP_genes,down=DOWN_genes)
	IsDE=list(up=1,down=-1)
	for (direction in c("up","down")) {
	  # zar : Διαφορικά εκφρασμένα ρυθμιστικά γονίδια
	  # (δηλ. που βρίσκονται στη στήλη 2 των αρχείων reference)
	  zar = intersect(TF[ ,2], Genes[[direction]])
    for (z in zar) { # z : parent (=1ου επιπέδου) γονίδιο
      # Κάθε εγγραφή περιέχει διαφορετικό πλήθος ρυθμιζόμενων γονιδίων,
      # που καθορίζεται στη στήλη 1 των αρχείων reference
      # και βρίσκονται από την στήλη 3 και πέρα
      last_index = 2 + TF[z, 1]
      regulated_genes = as.matrix(TF[z, 3:last_index])
      deup = intersect(regulated_genes, UP_genes)
      dedown = intersect(regulated_genes, DOWN_genes)
      deup_len = length(deup)
      dedown_len = length(dedown)
      if (needRegulatory) {
        # Parent - children
        if (deup_len > 0)
          Network = rbind(Network, cbind(z, deup))
        if (dedown_len > 0)
          Network = rbind(Network, cbind(z, dedown))
      }
      TF_counts[z, 3] = IsDE[[direction]]
      TF_counts[z, 4] = TF_counts[z, 4] + deup_len + dedown_len
      TF_counts[z, 5] = TF_counts[z, 5] + deup_len
      TF_counts[z, 6] = TF_counts[z, 6] + dedown_len

      #Για κάθε ρυθμιζόμενο γονίδιο
      for (tidx in 3:last_index) {
        gene = TF[z, tidx]
        if (length(intersect(gene, rownames(TF))) <= 0) next

        # αν ρυθμίζει κι αυτό κάποια γονίδια (είναι δηλαδή TF)
        inner_last_index=2 + TF[gene, 1]
        deup = intersect(as.matrix(TF[gene, 3:inner_last_index]), UP_genes)
        dedown = intersect(as.matrix(TF[gene, 3:inner_last_index]), DOWN_genes)
        deup_len = length(deup)
        dedown_len = length(dedown)
        if (needRegulatory) {
          # Parent - Child & Child - GrandChildren
          # Αν το ρυθμιζόμενο γονίδιο είναι DE, ίσως έχει ήδη προστεθεί
          # στο προηγούμενο στάδιο, αλλά είναι πιο εύκολο να αφαιρεθεί
          # η διπλή εγγραφή με unique() στο τέλος, από το να γίνονται
          # έλεγχοι εδώ
          if (deup_len > 0) {
            Network = rbind(Network, c(z, gene))
            Network = rbind(Network, cbind(gene, deup))
          }
          if (dedown_len > 0) {
            Network = rbind(Network, c(z, gene))
            Network = rbind(Network, cbind(gene, dedown))
          }
        }
        TF_counts[z, 1] = TF_counts[z, 1] + TF[gene, 1]
        TF_counts[z, 4] = TF_counts[z, 4] + deup_len + dedown_len
        TF_counts[z, 5] = TF_counts[z, 5] + deup_len
        TF_counts[z, 6] = TF_counts[z, 6] + dedown_len
      } # for (tidx...
    } # for (z...
	} # for (direction...

	for(mir in 1:nrow(miRNA)){
	  z=rownames(miRNA)[mir];
	  last_index=2+miRNA[z,1];
	  regulated_genes=as.matrix(miRNA[z,3:last_index])
	  deup=intersect(regulated_genes,UP_genes);
	  dedown=intersect(regulated_genes,DOWN_genes);
	  de=intersect(regulated_genes,DE_genes);
	  if(needRegulatory && length(de)>0)
	        Network=rbind(Network,cbind(z,de));
    # Πάει μόνο σε ένα επίπεδο, γιατί αν κάποιο από τα διαφορικά εκφρασμένα
	  # ρυθμιζόμενα γονίδια ήταν ρυθμιστής, θα είχαμε ήδη ασχοληθεί μαζί του
	  # στο προηγούμενο στάδιο
	  miRNA_counts[z,3]=miRNA_counts[z,3]+length(de);
	  miRNA_counts[z,4]=miRNA_counts[z,4]+length(deup);
	  miRNA_counts[z,5]=miRNA_counts[z,5]+length(dedown);
	} # for (mir...

	{
	# for(Tf in 1:nrow(TF)){
	# 	z<-intersect(TF[Tf,2],UP_genes);
	# 	if(length(z)>0){
	# 		deup<-intersect(as.matrix(TF[z,3:(2+TF[z,1])]),as.list(UP_genes));
	# 		dedown<-intersect(as.matrix(TF[z,3:(2+TF[z,1])]),as.list(DOWN_genes));
	# 		if(needRegulatory){
	# 			if(length(deup)>0){
	# 			  for(f in 1:length(deup)){
	# 					Network=rbind(Network,c(z,deup[f]));
	# 				}
	# 			}
	# 			if(length(dedown)>0){
	# 				for(f in 1:length(dedown)){
	# 					Network=rbind(Network,c(z,dedown[f]));
	# 				}
	# 			}
	# 		}
	# 		TF_counts[z,3]=1;
	# 		TF_counts[z,4]=TF_counts[z,4]+length(deup)+length(dedown);
	# 		TF_counts[z,5]=TF_counts[z,5]+length(deup);
	# 		TF_counts[z,6]=TF_counts[z,6]+length(dedown);
	# 		for(tar in 3:(TF[z,1]+2)){
	# 			if(length(intersect(TF[z,tar],rownames(TF)))>0){
	# 				deup<-intersect(as.matrix(TF[as.character(TF[z,tar]),3:(2+TF[as.character(TF[z,tar]),1])]),UP_genes);
	# 				dedown<-intersect(as.matrix(TF[as.character(TF[z,tar]),3:(2+TF[as.character(TF[z,tar]),1])]),DOWN_genes);
	# 				if(length(deup)>0){
	# 					for(f in 1:length(deup)){
	# 					  Network=rbind(Network,c(z,as.character(TF[z,tar])));
	# 						Network=rbind(Network,c(as.character(TF[z,tar]),deup[f]));
	# 					}
	# 				}
	# 				if(length(dedown)>0){
	# 					for(f in 1:length(dedown)){
	# 					  Network=rbind(Network,c(z,as.character(TF[z,tar])));
	# 						Network=rbind(Network,c(as.character(TF[z,tar]),dedown[f]));
	# 					}
	# 				}
	# 				TF_counts[z,1]=TF_counts[z,1]+TF[as.character(TF[z,tar]),1];
	# 				TF_counts[z,4]=TF_counts[z,4]+length(deup)+length(dedown);
	# 				TF_counts[z,5]=TF_counts[z,5]+length(deup);
	# 				TF_counts[z,6]=TF_counts[z,6]+length(dedown);
	# 			}
	# 		}
	# 	}
	# }
	#
	# for(Tf in 1:nrow(TF)){
	# 	z<-intersect(TF[Tf,2],DOWN_genes);
	# 	if(length(z)>0){
	# 	  deup<-intersect(as.matrix(TF[z,3:(2+TF[z,1])]),as.list(UP_genes));
	# 		dedown<-intersect(as.matrix(TF[z,3:(2+TF[z,1])]),as.list(DOWN_genes));
	# 		if((network=="global")||(network=="regulatory")){
	# 			if(length(deup)>0){
	# 				for(f in 1:length(deup)){
	# 					Network=rbind(Network,c(z,deup[f]));
	# 				}
	# 			}
	# 			if(length(dedown)>0){
	# 				for(f in 1:length(dedown)){
	# 					Network=rbind(Network,c(z,dedown[f]));
	# 				}
	# 			}
	# 		}
	# 		TF_counts[z,3]=-1;
	# 		TF_counts[z,4]=TF_counts[z,4]+length(deup)+length(dedown);
	# 		TF_counts[z,5]=TF_counts[z,5]+length(deup);
	# 		TF_counts[z,6]=TF_counts[z,6]+length(dedown);
	# 		for(tar in 3:(TF[z,1]+2)){
	# 			if(length(intersect(TF[z,tar],rownames(TF)))>0){
	# 				deup<-intersect(as.matrix(TF[as.character(TF[z,tar]),3:(2+TF[as.character(TF[z,tar]),1])]),UP_genes);
	# 				dedown<-intersect(as.matrix(TF[as.character(TF[z,tar]),3:(2+TF[as.character(TF[z,tar]),1])]),DOWN_genes);
	# 				if(length(deup)>0){
	# 					for(f in 1:length(deup)){
	# 						Network=rbind(Network,c(z,as.character(TF[z,tar])));
	# 						Network=rbind(Network,c(as.character(TF[z,tar]),deup[f]));
	# 					}
	# 				}
	# 				if(length(dedown)>0){
	# 					for(f in 1:length(dedown)){
	# 						Network=rbind(Network,c(z,as.character(TF[z,tar])));
	# 						Network=rbind(Network,c(as.character(TF[z,tar]),dedown[f]));
	# 					}
	# 				}
	# 				TF_counts[z,1]=TF_counts[z,1]+TF[as.character(TF[z,tar]),1];
	# 				TF_counts[z,4]=TF_counts[z,4]+length(deup)+length(dedown);
	# 				TF_counts[z,5]=TF_counts[z,5]+length(deup);
	# 				TF_counts[z,6]=TF_counts[z,6]+length(dedown);
	# 			}
	# 		}
	# 	}
	# }
  #
	# for(mir in 1:nrow(miRNA)){
	# 	z=rownames(miRNA)[mir];
	# 	deup<-intersect(as.matrix(miRNA[z,3:(2+miRNA[z,1])]),as.list(UP_genes));
	# 	dedown<-intersect(as.matrix(miRNA[z,3:(2+miRNA[z,1])]),as.list(DOWN_genes));
	# 	de<-intersect(as.matrix(miRNA[z,3:(2+miRNA[z,1])]),as.list(DE_genes));
	# 		if(needRegulatory){
	# 			if(length(de)>0){
	# 				for(f in 1:length(de)){
	# 					Network=rbind(Network,c(z,de[f]));
	# 				}
	# 			}
	# 		}
	# 	miRNA_counts[z,3]=miRNA_counts[z,3]+length(de);
	# 	miRNA_counts[z,4]=miRNA_counts[z,4]+length(deup);
	# 	miRNA_counts[z,5]=miRNA_counts[z,5]+length(dedown);
	# }

	# Κώστας (end)
	}

	for(keg in 1:nrow(kegg)){
		z=rownames(kegg)[keg];
		deup<-intersect(as.matrix(kegg[z,3:(2+kegg[z,1])]),as.list(UP_genes));
		dedown<-intersect(as.matrix(kegg[z,3:(2+kegg[z,1])]),as.list(DOWN_genes));
		de<-intersect(as.matrix(kegg[z,3:(2+kegg[z,1])]),as.list(DE_genes));
			if(needFunctional){
				if(length(de)>0){
					for(f in 1:length(de)){
						Network=rbind(Network,c(z,de[f]));
					}
				}
			}
		kegg_counts[z,3]=kegg_counts[z,3]+length(de);
		kegg_counts[z,4]=kegg_counts[z,4]+length(deup);
		kegg_counts[z,5]=kegg_counts[z,5]+length(dedown);
	}

	for(kegcat in 1:nrow(keggcat)){
		z=rownames(keggcat)[kegcat];
		deup<-intersect(as.matrix(keggcat[z,3:(2+keggcat[z,1])]),as.list(UP_genes));
		dedown<-intersect(as.matrix(keggcat[z,3:(2+keggcat[z,1])]),as.list(DOWN_genes));
		de<-intersect(as.matrix(keggcat[z,3:(2+keggcat[z,1])]),as.list(DE_genes));
			if(needFunctional){
				if(length(de)>0){
					for(f in 1:length(de)){
						Network=rbind(Network,c(z,de[f]));
					}
				}
			}
		keggcat_counts[z,3]=keggcat_counts[z,3]+length(de);
		keggcat_counts[z,4]=keggcat_counts[z,4]+length(deup);
		keggcat_counts[z,5]=keggcat_counts[z,5]+length(dedown);
	}

	for(go in 1:nrow(gos)){
		z=rownames(gos)[go];
		deup<-intersect(as.matrix(gos[z,3:(2+gos[z,1])]),as.list(UP_genes));
		dedown<-intersect(as.matrix(gos[z,3:(2+gos[z,1])]),as.list(DOWN_genes));
		de<-intersect(as.matrix(gos[z,3:(2+gos[z,1])]),as.list(DE_genes));
			if(needFunctional){
				if(length(de)>0){
					for(f in 1:length(de)){
						Network=rbind(Network,c(z,de[f]));
					}
				}
			}
		go_counts[z,3]=go_counts[z,3]+length(de);
		go_counts[z,4]=go_counts[z,4]+length(deup);
		go_counts[z,5]=go_counts[z,5]+length(dedown);
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
	    # Hypergeometric Distr.   # White_Drawn   Total_White                 Total_Black                  Total_Drawn
	                              # TargetsDE     TotalDE                     TotalNDE                     Targets
	    ResultsTFtemp[i,1]<-stats::phyper(TF_counts[r,4],Number_of_DE_genes,(TF_counts[r,2]-Number_of_DE_genes), TF_counts[r,1],lower.tail=FALSE);
	    ResultsTFtemp[i,2]<-stats::phyper(TF_counts[r,5],Number_of_Up_genes,(TF_counts[r,2]-Number_of_Up_genes), TF_counts[r,1],lower.tail=FALSE);
	    ResultsTFtemp[i,3]<-stats::phyper(TF_counts[r,6],Number_of_Down_genes,(TF_counts[r,2]-Number_of_Down_genes), TF_counts[r,1],lower.tail=FALSE);
	  }
	  ResultsTF[,1]=stats::p.adjust(ResultsTFtemp[,1],method="fdr");
	  ResultsTF[,2]=stats::p.adjust(ResultsTFtemp[,2],method="fdr");
	  ResultsTF[,3]=stats::p.adjust(ResultsTFtemp[,3],method="fdr");
	}

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

	if(length(which(kegg_counts[,3]!=0))==0){
		print("No Kegg pathway's member found deregulated..Really? Something has gone wrong.");
	} else {
# 		print(kegg_counts[which(kegg_counts[,3]!=0),]);
		Resultskegg=matrix(nrow=nrow(kegg_counts[which(kegg_counts[,3]!=0),]),ncol=3);
		rownames(Resultskegg)=rownames(kegg_counts)[which(kegg_counts[,3]!=0)];
		colnames(Resultskegg)=c("DE_qvalue","UP_qvalue","DOWN_qvalue");
		Resultskeggtemp=matrix(nrow=nrow(kegg_counts[which(kegg_counts[,3]!=0),]),ncol=3);
		rownames(Resultskeggtemp)=rownames(kegg_counts)[which(kegg_counts[,3]!=0)];
		colnames(Resultskeggtemp)=c("DE_pvalue","UP_pvalue","DOWN_pvalue");
		for(i in 1:length(which(kegg_counts[,3]!=0))){
			r=which(kegg_counts[,3]!=0)[i];
			Resultskeggtemp[i,1]<-stats::phyper(kegg_counts[r,3],Number_of_DE_genes,(kegg_counts[r,2]-Number_of_DE_genes), kegg_counts[r,1],lower.tail=FALSE);
			Resultskeggtemp[i,2]<-stats::phyper(kegg_counts[r,4],Number_of_Up_genes,(kegg_counts[r,2]-Number_of_Up_genes), kegg_counts[r,1],lower.tail=FALSE);
			Resultskeggtemp[i,3]<-stats::phyper(kegg_counts[r,5],Number_of_Down_genes,(kegg_counts[r,2]-Number_of_Down_genes), kegg_counts[r,1],lower.tail=FALSE);
		}
		Resultskegg[,1]=stats::p.adjust(Resultskeggtemp[,1],method="fdr");
		Resultskegg[,2]=stats::p.adjust(Resultskeggtemp[,2],method="fdr");
		Resultskegg[,3]=stats::p.adjust(Resultskeggtemp[,3],method="fdr");
	}

	if(length(which(keggcat_counts[,3]!=0))==0){
		print("No Kegg pathway category's member found deregulated..Really? Something has gone wrong.");
	} else {
# 		print(keggcat_counts[which(keggcat_counts[,3]!=0),]);
		Resultskeggcat=matrix(nrow=nrow(keggcat_counts[which(keggcat_counts[,3]!=0),]),ncol=3);
		rownames(Resultskeggcat)=rownames(keggcat_counts)[which(keggcat_counts[,3]!=0)];
		colnames(Resultskeggcat)=c("DE_qvalue","UP_qvalue","DOWN_qvalue");
		Resultskeggcattemp=matrix(nrow=nrow(keggcat_counts[which(keggcat_counts[,3]!=0),]),ncol=3);
		rownames(Resultskeggcattemp)=rownames(keggcat_counts)[which(keggcat_counts[,3]!=0)];
		colnames(Resultskeggcattemp)=c("DE_pvalue","UP_pvalue","DOWN_pvalue");
		for(i in 1:length(which(keggcat_counts[,3]!=0))){
			r=which(keggcat_counts[,3]!=0)[i];
			Resultskeggcattemp[i,1]<-stats::phyper(keggcat_counts[r,3],Number_of_DE_genes,(keggcat_counts[r,2]-Number_of_DE_genes), keggcat_counts[r,1],lower.tail=FALSE);
			Resultskeggcattemp[i,2]<-stats::phyper(keggcat_counts[r,4],Number_of_Up_genes,(keggcat_counts[r,2]-Number_of_Up_genes), keggcat_counts[r,1],lower.tail=FALSE);
			Resultskeggcattemp[i,3]<-stats::phyper(keggcat_counts[r,5],Number_of_Down_genes,(keggcat_counts[r,2]-Number_of_Down_genes), keggcat_counts[r,1],lower.tail=FALSE);
		}
		Resultskeggcat[,1]=stats::p.adjust(Resultskeggcattemp[,1],method="fdr");
		Resultskeggcat[,2]=stats::p.adjust(Resultskeggcattemp[,2],method="fdr");
		Resultskeggcat[,3]=stats::p.adjust(Resultskeggcattemp[,3],method="fdr");
	}

	if(length(which(go_counts[,3]!=0))==0){
		print("No GO's member found deregulated..Really? Something has gone wrong.");
	} else {
# 		print(head(go_counts[which(go_counts[,3]!=0),]));
		Resultsgo=matrix(nrow=nrow(go_counts[which(go_counts[,3]!=0),]),ncol=3);
		rownames(Resultsgo)=rownames(go_counts)[which(go_counts[,3]!=0)];
		colnames(Resultsgo)=c("DE_qvalue","UP_qvalue","DOWN_qvalue");
		Resultsgotemp=matrix(nrow=nrow(go_counts[which(go_counts[,3]!=0),]),ncol=3);
		rownames(Resultsgotemp)=rownames(go_counts)[which(go_counts[,3]!=0)];
		colnames(Resultsgotemp)=c("DE_pvalue","UP_pvalue","DOWN_pvalue");
		for(i in 1:length(which(go_counts[,3]!=0))){
			r=which(go_counts[,3]!=0)[i];
			Resultsgotemp[i,1]<-stats::phyper(go_counts[r,3],Number_of_DE_genes,(go_counts[r,2]-Number_of_DE_genes), go_counts[r,1],lower.tail=FALSE);
			Resultsgotemp[i,2]<-stats::phyper(go_counts[r,4],Number_of_Up_genes,(go_counts[r,2]-Number_of_Up_genes), go_counts[r,1],lower.tail=FALSE);
			Resultsgotemp[i,3]<-stats::phyper(go_counts[r,5],Number_of_Down_genes,(go_counts[r,2]-Number_of_Down_genes), go_counts[r,1],lower.tail=FALSE);
		}
		Resultsgo[,1]=stats::p.adjust(Resultsgotemp[,1],method="fdr");
		Resultsgo[,2]=stats::p.adjust(Resultsgotemp[,2],method="fdr");
		Resultsgo[,3]=stats::p.adjust(Resultsgotemp[,3],method="fdr");
	}

	# Network_final=matrix(ncol=2,data="");  # Κώστας - Όπως και πιο πάνω, δημιουργούσε μια κενή γραμμή
	Network_final=matrix(ncol=2,nrow=0);

	if(needFunctional){
	  # Για κάποιο λόγο από τον παλιό κώδικα η μεταβλητή Network είναι list,
	  # οπότε πρέπει να τη διορθώσω
	  Network=matrix(unlist(Network),ncol=2)
		for(interaction in 1:nrow(Network)){
			z1=intersect(Network[interaction,1],rownames(Resultskegg));
			z2=intersect(Network[interaction,1],rownames(Resultskeggcat));
			z3=intersect(Network[interaction,1],rownames(Resultsgo));
			if(length(z1>0)){
				if(Resultskegg[z1,1]<=PV_threshold){
					Network_final=rbind(Network_final,Network[interaction,]);
				}
			} else if(length(z2>0)){
				if(Resultskeggcat[z2,1]<=PV_threshold){
					Network_final=rbind(Network_final,Network[interaction,]);
				}
			} else if(length(z3>0)){
				if(Resultsgo[z3,1]<=PV_threshold){
					Network_final=rbind(Network_final,Network[interaction,]);
				}
			} else {
				Network_final=rbind(Network_final,Network[interaction,]);
			}
		}
	}
	if (needRegulatory) {
		Network_final=rbind(Network_final,Network);
	}
	Network_final=unique(Network_final)

	# Αν δεν υπάρχει το directory εξόδου, το δημιουργώ
	if (dir.exists(output_dir)==FALSE)
	  dir.create(output_dir)
	# Αν προϋπήρχε ή πέτυχε η δημιουργία το προσθέτω στο πρόθεμα
	if (dir.exists(output_dir)==TRUE)
	  output=file.path(output_dir,output)

	if (rank) {
	  P=RNRank(Network_final,...)
	  utils::write.csv(P,file=paste0(output,"_Ranks.csv"),quote = F,row.names = T)
	}
	utils::write.csv(Network_final,file=paste0(output,"_Network.csv"),quote=F,row.names=F);

	tempnames=c("Regulator",colnames(ResultsTF));
	ResultsTF=cbind(rownames(ResultsTF),as.data.frame(ResultsTF));
	colnames(ResultsTF)=tempnames;

	tempnames=c("miRNA",colnames(ResultsmiRNA));
	ResultsmiRNA=cbind(rownames(ResultsmiRNA),as.data.frame(ResultsmiRNA));
	colnames(ResultsmiRNA)=tempnames;

	tempnames=c("KEGG_pathway",colnames(Resultskegg));
	Resultskegg=cbind(rownames(Resultskegg),as.data.frame(Resultskegg));
	colnames(Resultskegg)=tempnames;

	tempnames=c("KEGG_pathway_category",colnames(Resultskeggcat));
	Resultskeggcat=cbind(rownames(Resultskeggcat),as.data.frame(Resultskeggcat));
	colnames(Resultskeggcat)=tempnames;

	tempnames=c("GO",colnames(Resultsgo));
	Resultsgo=cbind(rownames(Resultsgo),as.data.frame(Resultsgo));
	colnames(Resultsgo)=tempnames;


	if(type_of_output=="html"){
		sortable.html.table(as.data.frame(ResultsTF),paste(output,"TF_Enrichment.html"),page.title=paste(output,"TF_Enrichment"));
		sortable.html.table(as.data.frame(ResultsmiRNA),paste(output,"miRNA_Enrichment.html"),page.title=paste(output,"miRNA_Enrichment"));
		sortable.html.table(as.data.frame(Resultskegg),paste(output,"KEGG_Enrichment.html"),page.title=paste(output,"KEGG_Enrichment"));
		sortable.html.table(as.data.frame(Resultskeggcat),paste(output,"KEGG_category_Enrichment.html"),page.title=paste(output,"KEGG_category_Enrichment"));
		sortable.html.table(as.data.frame(Resultsgo),paste(output,"GO_Enrichment.html"),page.title=paste(output,"GO_Enrichment"));
	}else if(type_of_output=="csv"){
		utils::write.csv(as.data.frame(ResultsTF),file=paste0(output,"_TF_Enrichment.csv"),row.names=F);
		utils::write.csv(as.data.frame(ResultsmiRNA),file=paste0(output,"_miRNA_Enrichment.csv"),row.names=F);
		utils::write.csv(as.data.frame(Resultskegg),file=paste0(output,"_KEGG_Enrichment.csv"),row.names=F);
		utils::write.csv(as.data.frame(Resultskeggcat),file=paste0(output,"_KEGG_category_Enrichment.csv"),row.names=F);
		utils::write.csv(as.data.frame(Resultsgo),file=paste0(output,"_GO_Enrichment.csv"),row.names=F);
	}
	# Κώστας - Commented out
	# else {
	# 	print("Please select type of output to be either \"html\" or \"csv\"");
	# }

	print("Analysis completed!");
}



