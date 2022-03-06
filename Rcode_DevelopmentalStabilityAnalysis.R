
library("NSM3")
library("ppcor")
library("beeswarm")
library("biomaRt")
library("GSEABase")
library("GOstats")
library("AnnotationForge")
library("gplots")	
library("GO.db")
library("dplyr")
library("lattice")
library("maptools")
library("pvclust")
library("snow")
library("vioplot")

options(scipen=0)


PvalaueCutOff			<- 0.01
TPMCutOff				<- 0.1
TPMCutOff_TissueStage	<- log(2, base=10)

outputdir				<- paste("~/Documents/work/Research/Fluctuation/Data/Analysis/", Sys.Date(), "/", sep="")
dir.create(outputdir)
outputdir				<- paste(outputdir, "ALLinONE_LOG_RMcorrection/", sep="")
dir.create(outputdir)	
outputdir_tmp			<- paste(outputdir, "P-ValaueCutOff-", PvalaueCutOff, "_TPMCutOff-", TPMCutOff, "_TPMCutOffTissueStage-", TPMCutOff_TissueStage, "_NOTabsVARs/", sep = "")
dir.create(outputdir_tmp)




# Read gene expression tables (row=genes, column=samples consisting of embryos for 4 stages; TPM_all=for inbred, TPM_wild_all=for two wild populations, TPM_pool= for technical replicates)
	TPM_all				<- read.table("~/Documents/work/Research/Fluctuation/Data/ExpressionTables/TPM_NOmt_20M.txt", header=TRUE, sep="\t", row.names=1)
	gene_names			<- TPM_all[,1]
	names(gene_names)	<- rownames(TPM_all)
	name_inbred			<- colnames(TPM_all)

	TPM_all				<- TPM_all[,-1]
	TPM_all				<- log(TPM_all+1, base=10)
	allgene_list		<- rownames(TPM_all)

	TPM_wild_all		<- read.table("~/Documents/work/Research/Fluctuation/Data/ExpressionTables/TPM_wild_NOmt_20M.txt", header=TRUE, sep="\t", row.names=1)
	TPM_wild_all		<- TPM_wild_all[,-1]
	TPM_wild_all		<- log(TPM_wild_all+1, base=10)

	TPM_pool			<- read.table("~/Documents/work/Research/Fluctuation/Data/ExpressionTables/TPM_full_techreplicates_NOmt.txt", header=TRUE, sep="\t", row.names=1)
	TPM_pool			<- TPM_pool[,-1]
	TPM_pool			<- log(TPM_pool+1, base=10)




###################################################################################################################################
#### Functions for calculating developmental stability and evolutionary diversity (both for phenotype and gene expression levels) 
##################################################################################################################################

# For phenotypic level
if(1){
	# For developmental stability calculations. input=expression table (log-transformed TPM table) for each stage. output=returns an array of inter-individual distances.
	Make_dist_array			<- function(TPM_table){
		## Function to calculate the absolute difference in expression level of a gene. input= array of two individual expression levels of a gene.
		get_expdifference	<- function(TPM_pair){
			TPM_pair	<- unlist(TPM_pair)
			if(TPM_pair[1] == 0 && TPM_pair[2] == 0){
				return(NA)
			}else{
				if(TPM_pair[1] < TPMCutOff || TPM_pair[2] < TPMCutOff){
					return(NA)
				}else{
					return((TPM_pair[1] - TPM_pair[2]))
				}
			}
		}

		## In a given pair, use function "get_expdifference" for all genes and obtain the variance of the expression difference. -> for loop over all pairs.
		dist_array	<- numeric(ncol(TPM_table)/2)
		for (i in 1:(ncol(TPM_table)/2)){
			dist_array[i]	<- var(apply(TPM_table[,(2*i - 1):(2*i)], 1, get_expdifference), na.rm=TRUE)
		}
						
		return(dist_array)	
	}


	# For technical error calculations. Calculate technical error in the same method as for developmental stability
	Make_dist_array_tech	<- function(TPM_table){

		get_expdifference	<- function(TPM_pair){
			TPM_pair	<- unlist(TPM_pair)
			if(TPM_pair[1] == 0 && TPM_pair[2] == 0){
				return(NA)
			}else{
				if(TPM_pair[1] >= TPMCutOff && TPM_pair[2] >= TPMCutOff){
					return((TPM_pair[1] - TPM_pair[2]))
				}else{
					return(NA)
				}
			}
		}
		dist_array 		<- numeric(6)
		dist_array[1]	<- var(apply(TPM_table[,c(1,2)], 1, get_expdifference), na.rm=TRUE)
		dist_array[2]	<- var(apply(TPM_table[,c(1,3)], 1, get_expdifference), na.rm=TRUE)
		dist_array[3]	<- var(apply(TPM_table[,c(1,4)], 1, get_expdifference), na.rm=TRUE)
		dist_array[4]	<- var(apply(TPM_table[,c(2,3)], 1, get_expdifference), na.rm=TRUE)
		dist_array[5]	<- var(apply(TPM_table[,c(2,4)], 1, get_expdifference), na.rm=TRUE)
		dist_array[6]	<- var(apply(TPM_table[,c(3,4)], 1, get_expdifference), na.rm=TRUE)	

		return(dist_array)	
	}


	# For intraspecies diverisity calculations. Calculate intraspecific diversity in the same way as for developmental stability. input=expression table for wild population 1, expression table for wild population 2
	Make_dist_array_wild	<- function(TPM_table1, TPM_table2){
		
		get_expdifference	<- function(TPM_pair){
			TPM_pair	<- unlist(TPM_pair)
			if(TPM_pair[1] == 0 && TPM_pair[2] == 0){
				return(NA)
			}else{
				if(TPM_pair[1] >= TPMCutOff && TPM_pair[2] >= TPMCutOff){
					return((TPM_pair[1] - TPM_pair[2]))
				}else{
					return(NA)
				}
			}
		}
		dist_array	<- NA
		for (i in 1:ncol(TPM_table1)){
			for (j in 1:ncol(TPM_table2)){						
				dist_array	<- c(dist_array, var(apply(cbind(TPM_table1[,i], TPM_table2[,j]), 1, get_expdifference), na.rm=TRUE))
			}
		}	
		dist_array	<- dist_array[-1]
		
		return(dist_array)			
	}


	# Function to quantify developmental stability by 1-Spearman's rho (without threshold)
	Make_dist_array_spearman_0		<- function(TPM_table){
		num_pairs 	<- ncol(TPM_table) / 2
		dist_array	<- array()

		for (j in 1:num_pairs){
			if (i == 1){
				dist_array	<- 1 - cor(TPM_table[,(2*j - 1)], TPM_table[,2*j], method = "spearman")
			}else{
				tmp_dist 	<- 1 - cor(TPM_table[,(2*j - 1)], TPM_table[,2*j], method = "spearman")
				dist_array	<- c(dist_array, tmp_dist)
			}
		}
									
		return(dist_array)	
	}

	# Function to quantify developmental stability by 1-Spearman's rho (with threshold)
	Make_dist_array_spearman		<- function(TPM_table){
		num_pairs 	<- ncol(TPM_table) / 2
		dist_array	<- array()

		TPM_CutOff		<- function(TPM_value){
			if(TPM_value == 0){
				return (NA)
			}else{
				if(TPM_value >= TPMCutOff){
					return (TPM_value)
				}else{
					return (NA)
				}
			}
		}
		TPM_table	<- apply(TPM_table, c(1,2), TPM_CutOff)

		for (j in 1:num_pairs){
			if (i == 1){
				dist_array	<- 1 - cor(TPM_table[,(2*j - 1)], TPM_table[,2*j], method = "spearman", use="complete.obs")
			}else{
				tmp_dist 	<- 1 - cor(TPM_table[,(2*j - 1)], TPM_table[,2*j], method = "spearman", use="complete.obs")
				dist_array	<- c(dist_array, tmp_dist)
			}
		}
									
		return(dist_array)	
	}



	# Function to quantify intraspecies diversity by 1-Spearman's rho (without threshold)
	Make_dist_array_wild_spearman	<- function(TPM_table1, TPM_table2){
		dist_array	<- array()
		for (k in 1:ncol(TPM_table1)){
			for (j in 1:ncol(TPM_table2)){						
				if (k*j == 1){
					dist_array	<- 1 - cor(TPM_table1[,k], TPM_table2[,j], method = "spearman")
				}else{
					tmp_dist 	<- 1 - cor(TPM_table1[,k], TPM_table2[,j], method = "spearman")
					dist_array	<- c(dist_array, tmp_dist)
				}
			}
		}	

		return(dist_array)			
	}

}

# For gene-by-gene level
if(1){
	# Calculate stability for log-transformed gene expression arrays (no correction for gene expression level dependence). input=an array of expression levels for all inbred samples of a given developmental stage.
	get_fluctuation_log		<- function(TPM_table_row){
		TPM_table_row	<- unlist(TPM_table_row)
		
		array_flu <- numeric(length(TPM_table_row)/2)
		for (j in 1:(length(TPM_table_row)/2)){
			tmp_mean	<- mean(c(TPM_table_row[2*j-1], TPM_table_row[2*j]))		
			if(tmp_mean == 0){
				array_flu[j]	<- NA
			}else{
				if(TPM_table_row[2*j-1] >= TPMCutOff &&  TPM_table_row[2*j] >= TPMCutOff){
					array_flu[j]	<- abs(TPM_table_row[2*j-1] - TPM_table_row[2*j])
				}else{
					array_flu[j]	<- NA
				}
			}			
		}
		
		if(sum(is.na(array_flu)) == length(TPM_table_row)/2){
			fluctuation <- NA
		}else{
			fluctuation	<- mean(na.omit(array_flu))	
		}
	
		return(fluctuation)	
	}


	# Calculate technical error for log-transformed gene expression arrays (no correction for gene expression level dependence). input=an array of expression levels for all technical replicates of a given developmental stage.
	get_techerror_log		<- function(TPM_table_row){
		TPM_table_row	<- unlist(TPM_table_row)

		array_tech 		<- numeric(6)
		get_error_pair	<- function(embryo1, embryo2){
			tmp_mean	<- mean(c(embryo1, embryo2))
			if(tmp_mean == 0){
				return (NA)
			}else{
				if(embryo1 >= TPMCutOff &&  embryo2 >= TPMCutOff){
					return (abs(embryo1 - embryo2))
				}else{
					return (NA)
				}
			}			
		}
			
		array_tech[1]	<- get_error_pair(TPM_table_row[1], TPM_table_row[2])	
		array_tech[2]	<- get_error_pair(TPM_table_row[1], TPM_table_row[3])
		array_tech[3]	<- get_error_pair(TPM_table_row[1], TPM_table_row[4]) 
		array_tech[4]	<- get_error_pair(TPM_table_row[2], TPM_table_row[3])
		array_tech[5]	<- get_error_pair(TPM_table_row[2], TPM_table_row[4]) 
		array_tech[6]	<- get_error_pair(TPM_table_row[3], TPM_table_row[4])				
		
		if(sum(is.na(array_tech)) == 6){
			techerror <- NA
		}else{
			techerror<- mean(na.omit(array_tech))	
		}
	
		return(techerror)	
	}



	# Calculate technical error for log-transformed gene expression arrays (no correction for gene expression level dependence). input=an array of expression levels for all wild1 population samples + all wild2 population samples  of a given developmental stage.
	get_microevo_log2		<- function(TPM_table_row){
		wild1		<- unlist(TPM_table_row[1:(length(TPM_table_row)/2)])
		wild2		<- unlist(TPM_table_row[((length(TPM_table_row)/2)+1):length(TPM_table_row)])

		array_evo 	<- array()
		for (j in 1:length(wild1)){
			for (k in 1:length(wild2)){
				tmp_mean	<- mean(c(wild1[j], wild2[k]))		
				if(tmp_mean == 0){
					array_evo	<- c(array_evo, NA)
				}else{
					if(wild1[j] >= TPMCutOff &&  wild2[k] >= TPMCutOff){
						array_evo	<- c(array_evo, abs(wild1[j] - wild2[k]))
					}else{
						array_evo	<- c(array_evo, NA)
					}
				}			
			}
		}
		array_evo 	<- array_evo[-1]
		
		if(sum(is.na(array_evo)) == (length(wild1)*length(wild2))){
			microevo	<- NA
		}else{
			microevo	<- mean(na.omit(array_evo))	
		}
		
		return(microevo)	
	}


	# Calculate the interspecies diversity of expression levels of a gene. input=expression level of all samples of st.23.5 in medaka + array of expression levels of all samples of other species (at the phylotypic stage) 
	get_GXPchange_log		<- function(TPM_table_row){
		medaka	<- unlist(TPM_table_row[1:ncol(TPM_mid)])
		medaka	<- mean(medaka)
		other	<- unlist(TPM_table_row[(ncol(TPM_mid)+1):length(TPM_table_row)])
		other	<- mean(other)

		tmp_mean <- mean(c(medaka,other))
					
		if(tmp_mean == 0){
				change	<- NA
		}else{
				change	<- abs(medaka - other)
		}
		return(change)
	}

}


##############################################################################################################
#### Extract gene expression table by stage
##############################################################################################################

if(1){
	TPM_early	<- TPM_all[,grep("stage15", colnames(TPM_all))]
	TPM_mid		<- TPM_all[,grep("stage23.5", colnames(TPM_all))]
	TPM_late1	<- TPM_all[,grep("stage28", colnames(TPM_all))]
	TPM_late2	<- TPM_all[,grep("Hatch", colnames(TPM_all))]
	rm(TPM_all)
	gc()
	gc()


	TPM_early_pool			<- TPM_pool[,grep("stage15", colnames(TPM_pool))]
	TPM_mid_pool			<- TPM_pool[,grep("stage23.5", colnames(TPM_pool))]
	TPM_late1_pool			<- TPM_pool[,grep("stage28", colnames(TPM_pool))]
	TPM_late2_pool			<- TPM_pool[,grep("Hatch", colnames(TPM_pool))]
	rm(TPM_pool)
	gc()
	gc()


	TPM_Kasasa_early	<- TPM_wild_all[,grep("Ol_stage15_Kasasa", colnames(TPM_wild_all))]
	TPM_Kasasa_mid		<- TPM_wild_all[,grep("Ol_stage23.5_Kasasa", colnames(TPM_wild_all))]
	TPM_Kasasa_late1	<- TPM_wild_all[,grep("Ol_stage28_Kasasa", colnames(TPM_wild_all))]
	TPM_Kasasa_late2	<- TPM_wild_all[,grep("Ol_Hatch_Kasasa", colnames(TPM_wild_all))]		
	TPM_Oura_early		<- TPM_wild_all[,grep("Ol_stage15_Oura", colnames(TPM_wild_all))]
	TPM_Oura_mid		<- TPM_wild_all[,grep("Ol_stage23.5_Oura", colnames(TPM_wild_all))]
	TPM_Oura_late1		<- TPM_wild_all[,grep("Ol_stage28_Oura", colnames(TPM_wild_all))]
	TPM_Oura_late2		<- TPM_wild_all[,grep("Ol_Hatch_Oura", colnames(TPM_wild_all))]	
	rm(TPM_wild_all)
	gc()
	gc()
}



##############################################################################################################
#### Visualization of genes adopted for analysis
##############################################################################################################

if(1){
	outputdir_tmp2	<- paste(outputdir_tmp, "p-values/", sep="")
	dir.create(outputdir_tmp2)

	# Return colors for plots according to the p-value of test of difference between gene expression stability and technical error. output= return_flu (value of developmental stability), return_tech (value of technical error), return_color (color for plot), p-value
	get_color_array	<- function(TPMtable_row_bind){
			
		Fluctuation_row	<- unlist(TPMtable_row_bind[1:(length(TPMtable_row_bind)-4)])		
		Tech_row		<- unlist(TPMtable_row_bind[(length(TPMtable_row_bind)-3):length(TPMtable_row_bind)])		

		
		### Calculate tech error for each pair
		array_tech 		<- numeric(6)
		get_error_pair	<- function(embryo1, embryo2){
			tmp_mean	<- mean(c(embryo1, embryo2))
			if(tmp_mean == 0){
				return (NA)
			}else{
				return (abs(embryo1 - embryo2))
			}			
		}
			
		array_tech[1]	<- get_error_pair(Tech_row[1],Tech_row[2])
		array_tech[2]	<- get_error_pair(Tech_row[1],Tech_row[3])
		array_tech[3]	<- get_error_pair(Tech_row[1],Tech_row[4]) 
		array_tech[4]	<- get_error_pair(Tech_row[2],Tech_row[3])
		array_tech[5]	<- get_error_pair(Tech_row[2],Tech_row[4])
		array_tech[6]	<- get_error_pair(Tech_row[3],Tech_row[4])	

			
		### Calculate developmental stability for each pair
		array_flu <- numeric(length(Fluctuation_row)/2)   		
		for (j in 1:(length(Fluctuation_row)/2)){			
			tmp_mean	<- mean(c(Fluctuation_row[2*j-1], Fluctuation_row[2*j]))			
			if(tmp_mean == 0){
				array_flu[j]	<- NA
			}else{
				array_flu[j]	<-  abs(Fluctuation_row[2*j-1] - Fluctuation_row[2*j])	
			}
		}
			

		###　Test whether the mean of developmental stability is significantly greater than the mean of technical error			
		if(sum(is.na(array_tech)) == 6 || sum(is.na(array_flu)) == length(array_flu)){
			return_tech		<- NA
			return_flu		<- NA
			return_color	<- NA
			return_pvalue	<- NA						
		}else{
			tmp_ttest_result <- wilcox.test(array_flu, array_tech, alternative = "greater")
			if(tmp_ttest_result$p.value > PvalaueCutOff){
				return_tech		<- mean(na.omit(array_tech))
				return_flu		<- mean(na.omit(array_flu))
				return_color	<- rgb(127/255, 127/255, 127/255, alpha=0.3)
				return_pvalue	<- tmp_ttest_result$p.value	
			}else{
				return_tech		<- mean(na.omit(array_tech))
				return_flu		<- mean(na.omit(array_flu))
				return_color	<- "dodgerblue"
				return_pvalue	<- tmp_ttest_result$p.value					
			}	
		}

		return(c(return_flu, return_tech, return_color, return_pvalue))		
	}
	
	stages	<- c("early", "mid", "late1", "late2")
	for(i in 1:4){
		stage			<- stages[i]
		tmp_inbred		<- na.omit(get(paste("TPM_", stage, sep="")))
		tmp_tech		<- get(paste("TPM_", stage, "_pool", sep=""))[rownames(tmp_inbred),]
	
		tmp					<- na.omit(cbind(tmp_inbred, tmp_tech))
		flutechcolor_table	<- na.omit(t(data.frame(apply(tmp, 1, get_color_array))))


		## Visualization (developmental stability and technical error)
		if(1){
			tmp_fluctuation		<- flutechcolor_table[,1]
			tmp_techerror		<- flutechcolor_table[,2]
			color				<- flutechcolor_table[,3]
			
			tmp_input 			<- cbind(tmp_fluctuation[intersect(names(tmp_fluctuation), names(color))], tmp_techerror[intersect(names(tmp_fluctuation), names(color))])

			outputpath	<- paste(outputdir_tmp2, "fluctuation-techerror_", stage, ".pdf", sep="")
			pdf(file = outputpath, width = 7,	height = 7)
			par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 

			plot(	as.numeric(tmp_input[,1]), as.numeric(tmp_input[,2]), 
					xlim	= c(0.0001,2), 
					ylim	= c(0.0001,2), 
					#xaxt	= "n",
					#yaxt	= "n",
					log		= "xy",
					xlab	= "Fluctuation (log(E1)-log(E2) -> mean across pairs, no RM correlation)",
					ylab	= "Tech error (log(Ei)-log(Ek) -> mean across pairs, no RM correlation)",
					pch 	= 20,
					cex		= 0.3,
					col		= color,
					cex.lab = 0.8,
					main	= paste("Before gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
					las		= 1 
				)

			par(new = TRUE)

			yxfun	<- function(x){x}
			plot(	yxfun, 
					0.002, 2,
					xaxt="n", yaxt="n",
					xlab="", ylab=""　
				)

			genesubset	<- nrow(subset(flutechcolor_table, flutechcolor_table[,3]=="dodgerblue"))
			mtext(paste(genesubset, " genes", sep=""), side = 1, line = 6, adj = 1)
				
			dev.off()
		}


		## Visualization (developmental stability and expression level)
		if(1){
			tmp_fluctuation		<- flutechcolor_table[,1]
			tmp_explevel		<- tmp_inbred[rownames(tmp_fluctuation),]
			tmp_explevel		<- apply(tmp_inbred, 1, mean)
			
			tmp_input 			<- cbind(tmp_fluctuation[intersect(names(tmp_fluctuation), names(color))], tmp_explevel[intersect(names(tmp_fluctuation), names(color))])


			outputpath	<- paste(outputdir_tmp2, "fluctuation-explevel_", stage, ".pdf", sep="")
			pdf(file = outputpath, width = 7,	height = 7)
			par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 

			plot(	as.numeric(tmp_input[,1]), as.numeric(tmp_input[,2]), 
					xlim	= c(0.0001,5), 
					ylim	= c(0.0001,5), 
					#xaxt	= "n",
					#yaxt	= "n",
					log		= "xy",
					xlab	= "Fluctuation (log(E1)-log(E2) -> mean across pairs, no RM correlation)",
					ylab	= "Expression level (log -> mean across individuals)",
					pch 	= 20,
					cex		= 0.3,
					col		= color,
					cex.lab = 0.8,
					main	= paste("Before gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
					las		= 1 
				)

			par(new = TRUE)
			yxfun	<- function(x){x}
			plot(	yxfun, 
					0.002, 2,
					xaxt="n", yaxt="n",
					xlab="", ylab=""　
				)

			genesubset	<- nrow(subset(flutechcolor_table, flutechcolor_table[,3]=="dodgerblue"))
			mtext(paste(genesubset, " genes", sep=""), side = 1, line = 6, adj = 1)
				
			dev.off()
		}

		}

	rm(flutechcolor_table)
	gc()
	gc()

}



##############################################################################################################
#### Selection of genes for analysis
##############################################################################################################

if(1){
	genes_only_largeflu_PvalueCut	<- function(TPMtable_row_bind){
 	
		Fluctuation_row		<- unlist(TPMtable_row_bind[1:(length(TPMtable_row_bind)-4)])		
		Tech_row			<- unlist(TPMtable_row_bind[(length(TPMtable_row_bind)-3):length(TPMtable_row_bind)])				
					
		array_tech 		<- numeric(6)
		get_error_pair	<- function(embryo1, embryo2){
			tmp_mean	<- mean(c(embryo1, embryo2))
			if(tmp_mean == 0){
				return (NA)
			}else{
				return (abs(embryo1 - embryo2))
			}			
		}
			
		array_tech[1]	<- get_error_pair(Tech_row[1],Tech_row[2])
		array_tech[2]	<- get_error_pair(Tech_row[1],Tech_row[3])
		array_tech[3]	<- get_error_pair(Tech_row[1],Tech_row[4])
		array_tech[4]	<- get_error_pair(Tech_row[2],Tech_row[3])
		array_tech[5]	<- get_error_pair(Tech_row[2],Tech_row[4])
		array_tech[6]	<- get_error_pair(Tech_row[3],Tech_row[4])	

				
		array_flu <- numeric(length(Fluctuation_row)/2)   		
		for (j in 1:(length(Fluctuation_row)/2)){			
			tmp_mean	<- mean(c(Fluctuation_row[2*j-1], Fluctuation_row[2*j]))			
			if(tmp_mean == 0){
				array_flu[j]	<- NA
			}else{
				array_flu[j]	<- abs(Fluctuation_row[2*j-1] - Fluctuation_row[2*j])	
			}
		}
			
			
		if(sum(is.na(array_tech)) == 6 || sum(is.na(array_flu)) == length(array_flu)){
			return_genes	<- rep(NA, length(array_flu))				
		}else{
			tmp_ttest_result <- wilcox.test(array_flu, array_tech, alternative = "greater")
			if(tmp_ttest_result$p.value > PvalaueCutOff){
				return_genes	<- rep(NA, length(Fluctuation_row))
			}else{
				return_genes	<- Fluctuation_row		
			}
		}
		return(return_genes)
	}



	stages	<- c("early", "mid", "late1", "late2")
	for(i in 1:4){
		stage			<- stages[i]

		tmp_inbred		<- na.omit(get(paste("TPM_", stage, sep="")))	
		tmp_tech		<- get(paste("TPM_", stage, "_pool", sep=""))[rownames(tmp_inbred),]
		
		tmp				<- na.omit(cbind(tmp_inbred, tmp_tech)) 
		tmp				<- t(data.frame(apply(tmp, 1, genes_only_largeflu_PvalueCut)))
		assign(paste("TPM_", stage, sep=""), na.omit(tmp))
	}		
		
	TPM_Kasasa_early	<- TPM_Kasasa_early[rownames(na.omit(TPM_early)),]
	TPM_Kasasa_mid		<- TPM_Kasasa_mid[rownames(na.omit(TPM_mid)),]
	TPM_Kasasa_late1	<- TPM_Kasasa_late1[rownames(na.omit(TPM_late1)),]
	TPM_Kasasa_late2	<- TPM_Kasasa_late2[rownames(na.omit(TPM_late2)),]		

	TPM_Oura_early		<- TPM_Oura_early[rownames(na.omit(TPM_early)),]
	TPM_Oura_mid		<- TPM_Oura_mid[rownames(na.omit(TPM_mid)),]
	TPM_Oura_late1		<- TPM_Oura_late1[rownames(na.omit(TPM_late1)),]
	TPM_Oura_late2		<- TPM_Oura_late2[rownames(na.omit(TPM_late2)),]
}


##############################################################################################################
#### Quantification of developmental stability at the phenotypic level, and comparison between stages
##############################################################################################################

if(1){
	outputdir_tmp2	<- paste(outputdir_tmp, "stage-boxplots_AfterGeneSelection/", sep="")
	dir.create(outputdir_tmp2)
				
	# Developmental stability
	if(1){	
		dist_early	<- Make_dist_array(TPM_early)
		dist_mid	<- Make_dist_array(TPM_mid)
		dist_late1	<- Make_dist_array(TPM_late1)
		dist_late2	<- Make_dist_array(TPM_late2)

		outputpath	<- paste(outputdir_tmp2, "Box_inbred.pdf", sep="")
		pdf(file = outputpath,	width = 7, height = 7)
   		par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1,dist_late2,
									horizontal	= TRUE,
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					#xaxt	= "n",
									ylim	= c(0, 0.06),
				 					cex.axis= 1,
									cex 	= 0.3,
									pch 	= 16,
									main	= paste("After gene selection, Log, inbred, genes with TPM=>",TPMCutOff, sep=""),
									las		= 1
								)

		axis(	side		= 4,
				at			= 1:4,
    			labels		= c(paste( length(dist_early), " pairs\n(n = ", ncol(TPM_early), ")",sep=""),
    							paste( length(dist_mid), " pairs\n(n = ", ncol(TPM_mid), ")",sep=""),
    							paste( length(dist_late1), " pairs\n(n = ", ncol(TPM_late1), ")",sep=""),
    							paste( length(dist_late2), " pairs\n(n = ", ncol(TPM_late2), ")",sep="")
    							),
    			tck			= 0,
    			las			= 1,
   		 		mgp			= c(3.5,0.7,0),
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)

		par(new=T)
			
		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))

		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xlim		= c(0, 0.06),
					xaxt		= "n",
					yaxt		= "n",
					ylab		= "",
					xlab		= "Variance of TPM differential",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5)
									)   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
					)

		# Statistical test Kruskal Wallis)
		KW_result	<- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")	
		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()
	}


	# Developmental stability and technical error
	if(1){	
		outputpath	<- paste(outputdir_tmp2, "Box_InbredandTech.pdf", sep="")
		pdf(file = outputpath,	width = 7, height = 7)
   		par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		dist_early	<- Make_dist_array(TPM_early)
		dist_mid	<- Make_dist_array(TPM_mid)
		dist_late1	<- Make_dist_array(TPM_late1)
		dist_late2	<- Make_dist_array(TPM_late2)

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1,dist_late2,
									horizontal	= TRUE,
									col		= "deepskyblue",
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					ylim	= c(0, 0.06),
									xlab		= "Variance of TPM differential",
				 					cex.axis= 1,
									cex 	= 0.3,
									pch 	= 16,
									main	= paste("After gene selection, Log, inbred, genes with TPM=>",TPMCutOff, ", inbred", sep=""),
									las		= 1
								)

		par(new=T)

		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim		= c(0, 0.06),
					ylab		= "",
					xlab		= "",
					cex			= 0.9,
					pch 		= 16, 
					col 		= c(bg=rgb(0,191/255,255/255, alpha=0.5))
					)

		par(new=T)

		dist_early_tech	<- Make_dist_array_tech(TPM_early_pool)
		dist_mid_tech	<- Make_dist_array_tech(TPM_mid_pool)
		dist_late1_tech	<- Make_dist_array_tech(TPM_late1_pool)
		dist_late2_tech	<- Make_dist_array_tech(TPM_late2_pool)

		box_position	<- boxplot(	dist_early_tech, dist_mid_tech, dist_late1_tech,dist_late2_tech,
									horizontal	= TRUE,
									col		= "gray50",
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					ylim	= c(0, 0.06),
				 					xaxt	= "n",
									yaxt	= "n",
				 					cex.axis= 1,
									cex 	= 0.3,
									pch 	= 16,
									las		= 1
								)

		par(new=T)

		dist				<- c(dist_early_tech, dist_mid_tech, dist_late1_tech, dist_late2_tech)
		stage				<- c(rep("Early",length(dist_early_tech)),rep("Mid",length(dist_mid_tech)),rep("Late1",length(dist_late1_tech)),rep("Late2",length(dist_late2_tech)))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim		= c(0, 0.06),
					ylab		= "",
					xlab		= "",
					cex			= 0.9,
					pch 		= 16, 
					col 		= c(bg=rgb(127/255,127/255,127/255, alpha=0.5)) 
					)

		axis(	side		= 4,
				at			= 1:4,
    			labels		= c(paste("Median\n", "inbred:", signif(median(dist_early), digits=3), ", tech:", signif(median(dist_early_tech), digits=3), sep=""),
    							paste("Median\n", "inbred:", signif(median(dist_mid), digits=3), ", tech:", signif(median(dist_mid_tech), digits=3), sep=""),
    							paste("Median\n", "inbred:", signif(median(dist_late1), digits=3), ", tech:", signif(median(dist_late1_tech), digits=3), sep=""),
    							paste("Median\n", "inbred:", signif(median(dist_late2), digits=3), ", tech:", signif(median(dist_late2_tech), digits=3), sep="")
    							),
    			tck			= 0,
    			las			= 1, 
   		 		mgp			= c(3.5,0.7,0), 
   		 		cex.axis	= 0.65,
   		 		adj			= 0.5
			)

		# Statistical test (Kruskal Wallis)
		KW_result	<- kruskal.test(x=list(dist_early_tech,dist_mid_tech,dist_late1_tech,dist_late2_tech))
		mtext(paste("among tech replicates, P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex=0.9)

		# Statistical test (one-sided wilcox)
		tmp	<- wilcox.test(dist_early, dist_early_tech, alternative = "greater")
		mtext(paste("inbred vs tech, P = ", signif(tmp$p.val, 3), " (early, one-sided wilcox)", sep=""), side = 1, line = 7, adj = 1, cex=0.9)
		tmp	<- wilcox.test(dist_mid, dist_mid_tech, alternative = "greater")
		mtext(paste("inbred vs tech, P = ", signif(tmp$p.val, 3), " (mid, one-sided wilcox)", sep=""), side = 1, line = 8, adj = 1, cex=0.9)
		tmp	<- wilcox.test(dist_late1, dist_late1_tech, alternative = "greater")
		mtext(paste("inbred vs tech, P = ", signif(tmp$p.val, 3), " (late1, one-sided wilcox)", sep=""), side = 1, line = 9, adj = 1, cex=0.9)
		tmp	<- wilcox.test(dist_late2, dist_late2_tech, alternative = "greater")
		mtext(paste("inbred vs tech, P = ", signif(tmp$p.val, 3), " (late2, one-sided wilcox)", sep=""), side = 1, line = 10, adj = 1, cex=0.9)

		dev.off()
	}


	# Intraspecies diversity (wild Kasasa vs wild Oura)
	if(1){
		dist_early			<- Make_dist_array_wild(TPM_Kasasa_early, TPM_Oura_early)
		dist_mid			<- Make_dist_array_wild(TPM_Kasasa_mid, TPM_Oura_mid)
		dist_late1			<- Make_dist_array_wild(TPM_Kasasa_late1, TPM_Oura_late1)
		dist_late2			<- Make_dist_array_wild(TPM_Kasasa_late2, TPM_Oura_late2)

		outputpath	<- paste(outputdir_tmp2, "Box_wild.pdf", sep="")
		pdf(file = outputpath,	width = 7,		height = 7)
    	par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1, dist_late2,
									horizontal	= TRUE,
									ylim	= c(0, 0.1),
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					cex.axis	= 1,
									cex 	= 0.3,
									pch 	= 16,
									main 	= paste("After gene selection, Log, wild, genes with TPM=>",TPMCutOff, sep=""),
									las		= 1
								)

		par(new=T)

		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim	= c(0, 0.1),
					ylab		= "",
					xlab		= "Variance of TPM differential",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5))   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
				)

		axis(	side		= 4, 
				at			= 1:4,
    			labels		= c(paste( length(dist_early), " pairs\n(n = 4, 4)",sep=""),
    							paste( length(dist_mid), " pairs\n(n = 4, 4)",sep=""),
    							paste( length(dist_late1), " pairs\n(n = 4, 4)",sep=""),
    							paste( length(dist_late2), " pairs\n(n = 2, 2)",sep="")
    							),
    			tck			= 0,
    			las			= 1, 
   		 		mgp			= c(3.5,0.7,0),
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)
		
		# Statistical test Kruskal Wallis)
		KW_result <- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")

		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()

	}


	# Intraspecies diversity (inbred vs wild Kasasa)
	if(1){
		TPM_early_F	<- TPM_early[,c(1:8,23,24,37:40)]
		TPM_mid_F	<- TPM_mid[,c(1:8,15,16,19,20,25,26,33,34,37:42,45,46)]
		TPM_late1_F	<- TPM_late1[,c(1:10,17,18,21:26,31,32,37,38,43,44,47,48)]
		TPM_late2_F	<- TPM_late2[,c(3:16,21:26)]

		dist_early	<- Make_dist_array_wild(TPM_Kasasa_early, TPM_early_F)
		dist_mid	<- Make_dist_array_wild(TPM_Kasasa_mid, TPM_mid_F)
		dist_late1	<- Make_dist_array_wild(TPM_Kasasa_late1, TPM_late1_F)
		dist_late2	<- Make_dist_array_wild(TPM_Kasasa_late2, TPM_late2_F)

		outputpath	<- paste(outputdir_tmp2, "Box_inbred-Kasasa.pdf", sep="")
		pdf(file = outputpath,	width = 7,		height = 7)
    	par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1, dist_late2,
									horizontal	= TRUE,
									ylim	= c(0, 0.2),
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					cex.axis	= 1,
									cex 	= 0.3,
									pch 	= 16,
									main 	= paste("After gene selection, Log, inbred-Kasasa, genes with TPM=>",TPMCutOff, sep=""),
									las		= 1
								)

		par(new=T)

		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim	= c(0, 0.2),
					ylab		= "",
					xlab		= "Variance of TPM differential",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5))   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
				)

		axis(	side		= 4, 
				at			= 1:4, 
    			labels		= c(paste( length(dist_early), " pairs\n(n = ", ncol(TPM_Kasasa_early),  ", ", ncol(TPM_early_F), ")",sep=""),
    							paste( length(dist_mid), " pairs\n(n = ", ncol(TPM_Kasasa_mid),  ", ", ncol(TPM_mid_F), ")",sep=""),
    							paste( length(dist_late1), " pairs\n(n = ", ncol(TPM_Kasasa_late1),  ", ", ncol(TPM_late1_F), ")",sep=""),
    							paste( length(dist_late2), " pairs\n(n = ", ncol(TPM_Kasasa_late2),  ", ", ncol(TPM_late2_F), ")",sep="")
    							),
    			tck			= 0,
    			las			= 1,  
   		 		mgp			= c(3.5,0.7,0), 
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)
		
		# Statistical test Kruskal Wallis)
		KW_result <- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")

		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()

	}


	# Intraspecies diversity (inbred vs wild Oura)
	if(1){
		TPM_early_F	<- TPM_early[,c(1:8,23,24,37:40)]
		TPM_mid_F	<- TPM_mid[,c(1:8,15,16,19,20,25,26,33,34,37:42,45,46)]
		TPM_late1_F	<- TPM_late1[,c(1:10,17,18,21:26,31,32,37,38,43,44,47,48)]
		TPM_late2_F	<- TPM_late2[,c(3:16,21:26)]

		dist_early	<- Make_dist_array_wild(TPM_early_F, TPM_Oura_early)
		dist_mid	<- Make_dist_array_wild(TPM_mid_F, TPM_Oura_mid)
		dist_late1	<- Make_dist_array_wild(TPM_late1_F, TPM_Oura_late1)
		dist_late2	<- Make_dist_array_wild(TPM_late2_F, TPM_Oura_late2)

		outputpath	<- paste(outputdir_tmp2, "Box_inbred-Oura.pdf", sep="")
		pdf(file = outputpath,	width = 7,		height = 7)
    	par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1, dist_late2,
									horizontal	= TRUE,
									ylim	= c(0, 0.2),
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					cex.axis	= 1,
									cex 	= 0.3,
									pch 	= 16,
									main 	= paste("After gene selection, Log, inbred-Oura, genes with TPM=>",TPMCutOff, sep=""),
									las		= 1
								)

		par(new=T)

		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim	= c(0, 0.2),
					ylab		= "",
					xlab		= "Variance of TPM differential",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5))   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
				)

		axis(	side		= 4, 
				at			= 1:4,  
    			labels		= c(paste( length(dist_early), " pairs\n(n = ", ncol(TPM_Oura_early),  ", ", ncol(TPM_early_F), ")",sep=""),
    							paste( length(dist_mid), " pairs\n(n = ", ncol(TPM_Oura_mid),  ", ", ncol(TPM_mid_F), ")",sep=""),
    							paste( length(dist_late1), " pairs\n(n = ", ncol(TPM_Oura_late1),  ", ", ncol(TPM_late1_F), ")",sep=""),
    							paste( length(dist_late2), " pairs\n(n = ", ncol(TPM_Oura_late2),  ", ", ncol(TPM_late2_F), ")",sep="")
    							),
    			tck			= 0,
    			las			= 1,  	
   		 		mgp			= c(3.5,0.7,0), 	
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)
		
		# Statistical test Kruskal Wallis)
		KW_result <- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")

		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()

	}


	# Intraspecies diversity (inbred vs wild Kasasa, by 1- Spearman's rho)
	if(1){
		TPM_early_F	<- TPM_early[,c(1:8,23,24,37:40)]
		TPM_mid_F	<- TPM_mid[,c(1:8,15,16,19,20,25,26,33,34,37:42,45,46)]
		TPM_late1_F	<- TPM_late1[,c(1:10,17,18,21:26,31,32,37,38,43,44,47,48)]
		TPM_late2_F	<- TPM_late2[,c(3:16,21:26)]

		dist_early	<- Make_dist_array_wild_spearman(TPM_Kasasa_early, TPM_early_F)
		dist_mid	<- Make_dist_array_wild_spearman(TPM_Kasasa_mid, TPM_mid_F)
		dist_late1	<- Make_dist_array_wild_spearman(TPM_Kasasa_late1, TPM_late1_F)
		dist_late2	<- Make_dist_array_wild_spearman(TPM_Kasasa_late2, TPM_late2_F)

		outputpath	<- paste(outputdir_tmp2, "Box_inbred-Kasasa_spearman.pdf", sep="")
		pdf(file = outputpath,	width = 7,		height = 7)
    	par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1, dist_late2,
									horizontal	= TRUE,
									ylim	= c(0, 0.3),
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					cex.axis	= 1,
									cex 	= 0.3,
									pch 	= 16,
									main 	= paste("After gene selection, Log, inbred-Kasasa, genes with TPM=>",TPMCutOff, sep=""),
									las		= 1
								)

		par(new=T)

		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim	= c(0, 0.3),
					ylab		= "",
					xlab		= "1 - Spearman's rho",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5))   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
				)

		axis(	side		= 4, 
				at			= 1:4,  
    			labels		= c(paste( length(dist_early), " pairs\n(n = ", ncol(TPM_Kasasa_early),  ", ", ncol(TPM_early_F), ")",sep=""),
    							paste( length(dist_mid), " pairs\n(n = ", ncol(TPM_Kasasa_mid),  ", ", ncol(TPM_mid_F), ")",sep=""),
    							paste( length(dist_late1), " pairs\n(n = ", ncol(TPM_Kasasa_late1),  ", ", ncol(TPM_late1_F), ")",sep=""),
    							paste( length(dist_late2), " pairs\n(n = ", ncol(TPM_Kasasa_late2),  ", ", ncol(TPM_late2_F), ")",sep="")
    							),
    			tck			= 0,		
    			las			= 1,  			
   		 		mgp			= c(3.5,0.7,0), 	
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)
		
		# Statistical test Kruskal Wallis)
		KW_result <- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")

		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()

	}


	# Intraspecies diversity (inbred vs wild Oura, by 1- Spearman's rho)
	if(1){
		TPM_early_F	<- TPM_early[,c(1:8,23,24,37:40)]
		TPM_mid_F	<- TPM_mid[,c(1:8,15,16,19,20,25,26,33,34,37:42,45,46)]
		TPM_late1_F	<- TPM_late1[,c(1:10,17,18,21:26,31,32,37,38,43,44,47,48)]
		TPM_late2_F	<- TPM_late2[,c(3:16,21:26)]

		dist_early	<- Make_dist_array_wild_spearman(TPM_early_F, TPM_Oura_early)
		dist_mid	<- Make_dist_array_wild_spearman(TPM_mid_F, TPM_Oura_mid)
		dist_late1	<- Make_dist_array_wild_spearman(TPM_late1_F, TPM_Oura_late1)
		dist_late2	<- Make_dist_array_wild_spearman(TPM_late2_F, TPM_Oura_late2)

		outputpath	<- paste(outputdir_tmp2, "Box_inbred-Oura_spearman.pdf", sep="")
		pdf(file = outputpath,	width = 7,		height = 7)
    	par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1, dist_late2,
									horizontal	= TRUE,
									ylim	= c(0, 0.25),
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					cex.axis	= 1,
									cex 	= 0.3,
									pch 	= 16,
									main 	= paste("After gene selection, Log, inbred-Oura, genes with TPM=>",TPMCutOff, sep=""),
									las		= 1
								)

		par(new=T)

		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim	= c(0, 0.25),
					ylab		= "",
					xlab		= "1 - Spearman's rho",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5))   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
				)

		axis(	side		= 4, 		
				at			= 1:4,  	
    			labels		= c(paste( length(dist_early), " pairs\n(n = ", ncol(TPM_Oura_early),  ", ", ncol(TPM_early_F), ")",sep=""),
    							paste( length(dist_mid), " pairs\n(n = ", ncol(TPM_Oura_mid),  ", ", ncol(TPM_mid_F), ")",sep=""),
    							paste( length(dist_late1), " pairs\n(n = ", ncol(TPM_Oura_late1),  ", ", ncol(TPM_late1_F), ")",sep=""),
    							paste( length(dist_late2), " pairs\n(n = ", ncol(TPM_Oura_late2),  ", ", ncol(TPM_late2_F), ")",sep="")
    							),
    			tck			= 0,		
    			las			= 1,  			
   		 		mgp			= c(3.5,0.7,0), 
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)
		
		# Statistical test Kruskal Wallis)
		KW_result <- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")

		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()

	}


	# Developmental stability (only female samples)
	if(1){	
		TPM_early_F	<- TPM_early[,c(1:8,23,24,37:40)]
		TPM_mid_F	<- TPM_mid[,c(1:8,15,16,19,20,25,26,33,34,37:42,45,46)]
		TPM_late1_F	<- TPM_late1[,c(1:10,17,18,21:26,31,32,37,38,43,44,47,48)]
		TPM_late2_F	<- TPM_late2[,c(3:16,21:26)]

		dist_early	<- Make_dist_array(TPM_early_F)
		dist_mid	<- Make_dist_array(TPM_mid_F)
		dist_late1	<- Make_dist_array(TPM_late1_F)
		dist_late2	<- Make_dist_array(TPM_late2_F)

		outputpath	<- paste(outputdir_tmp2, "Box_inbredFemale.pdf", sep="")
		pdf(file = outputpath,	width = 7, height = 7)
   		par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1,dist_late2,
									horizontal	= TRUE,
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					ylim	= c(0, 0.06),
				 					#xaxt	= "n",
				 					cex.axis= 1,
									cex 	= 0.3,
									pch 	= 16,
									main	= paste("After gene selection, Log, inbred Female, genes with TPM=>",TPMCutOff, ", inbred", sep=""),
									las		= 1
								)

		axis(	side		= 4,
				at			= 1:4,  
    			labels		= c(paste( length(dist_early), " pairs\n(n = ", ncol(TPM_early_F), ")",sep=""),
    							paste( length(dist_mid), " pairs\n(n = ", ncol(TPM_mid_F), ")",sep=""),
    							paste( length(dist_late1), " pairs\n(n = ", ncol(TPM_late1_F), ")",sep=""),
    							paste( length(dist_late2), " pairs\n(n = ", ncol(TPM_late2_F), ")",sep="")
    							),
    			tck			= 0,	
    			las			= 1,  		
   		 		mgp			= c(3.5,0.7,0), 	
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)

		par(new=T)
			
		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))

		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim	= c(0, 0.06),
					ylab		= "",
					xlab		= "Variance of TPM differential",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5)
									)   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
					)

		# Statistical test Kruskal Wallis)
		KW_result	<- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")	
		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()
	}


	# Developmental stability (only male samples)
	if(1){	
		TPM_early_M	<- TPM_early[,c(1:8,23,24,37:40)*(-1)]
		TPM_mid_M	<- TPM_mid[,c(1:8,15,16,19,20,25,26,33,34,37:42,45,46)*(-1)]
		TPM_late1_M	<- TPM_late1[,c(1:10,17,18,21:26,31,32,37,38,43,44,47,48)*(-1)]
		TPM_late2_M	<- TPM_late2[,c(3:16,21:26)*(-1)]

		dist_early	<- Make_dist_array(TPM_early_M)
		dist_mid	<- Make_dist_array(TPM_mid_M)
		dist_late1	<- Make_dist_array(TPM_late1_M)
		dist_late2	<- Make_dist_array(TPM_late2_M)

		outputpath	<- paste(outputdir_tmp2, "Box_inbredMale.pdf", sep="")
		pdf(file = outputpath,	width = 7, height = 7)
   		par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1,dist_late2,
									horizontal	= TRUE,
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					ylim	= c(0, 0.06),
				 					#xaxt	= "n",
				 					cex.axis= 1,
									cex 	= 0.3,
									pch 	= 16,
									main	= paste("After gene selection, Log, inbred Male, genes with TPM=>",TPMCutOff, ", inbred", sep=""),
									las		= 1
								)

		axis(	side		= 4,
				at			= 1:4,  
    			labels		= c(paste( length(dist_early), " pairs\n(n = ", ncol(TPM_early_M), ")",sep=""),
    							paste( length(dist_mid), " pairs\n(n = ", ncol(TPM_mid_M), ")",sep=""),
    							paste( length(dist_late1), " pairs\n(n = ", ncol(TPM_late1_M), ")",sep=""),
    							paste( length(dist_late2), " pairs\n(n = ", ncol(TPM_late2_M), ")",sep="")
    							),
    			tck			= 0,	
    			las			= 1,  		
   		 		mgp			= c(3.5,0.7,0), 	
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)

		par(new=T)
			
		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))

		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim	= c(0, 0.06),
					ylab		= "",
					xlab		= "Variance of TPM differential",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5)
									)   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
					)

		# Statistical test Kruskal Wallis)
		KW_result	<- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")	
		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()
	}


	# Developmental stability (with genes shared between stages)
	if(0){	
		TPM_early_shared	<- na.omit(TPM_early)
		TPM_mid_shared		<- na.omit(TPM_mid)
		TPM_late1_shared	<- na.omit(TPM_late1)
		TPM_late2_shared	<- na.omit(TPM_late2)
		
		shared_genelist		<- intersect(rownames(TPM_early_shared), rownames(TPM_mid_shared))
		shared_genelist		<- intersect(shared_genelist, rownames(TPM_late1_shared))
		shared_genelist		<- intersect(shared_genelist, rownames(TPM_late2_shared))

		dist_early			<- Make_dist_array(TPM_early[shared_genelist,])
		dist_mid			<- Make_dist_array(TPM_mid[shared_genelist,])
		dist_late1			<- Make_dist_array(TPM_late1[shared_genelist,])
		dist_late2			<- Make_dist_array(TPM_late2[shared_genelist,])

		outputpath	<- paste(outputdir_tmp2, "Box_inbred_samegeneset.pdf", sep="")
		pdf(file = outputpath,	width = 7, height = 7)
   		par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1,dist_late2,
									horizontal	= TRUE,
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					ylim	= c(0, 0.15),
				 					#xaxt	= "n",
				 					cex.axis= 1,
									cex 	= 0.3,
									pch 	= 16,
									main	= paste("After gene selection, Log, inbred(only shared genes), genes with TPM=>",TPMCutOff, ", inbred", sep=""),
									las		= 1
								)

		axis(	side		= 4,
				at			= 1:4,  		# 座標
    			labels		= c(paste( length(dist_early), " pairs\n(n = ", ncol(TPM_early), ")",sep=""),
    							paste( length(dist_mid), " pairs\n(n = ", ncol(TPM_mid), ")",sep=""),
    							paste( length(dist_late1), " pairs\n(n = ", ncol(TPM_late1), ")",sep=""),
    							paste( length(dist_late2), " pairs\n(n = ", ncol(TPM_late2), ")",sep="")
    							),
    			tck			= 0,				#　目盛りの長さ、内側に食い込む 
    			las			= 1,  				# ラベルスタイル（las）は1（常に水平）
   		 		mgp			= c(3.5,0.7,0), 		# ラベル位置
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)

		par(new=T)
			
		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))

		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim	= c(0, 0.15),
					ylab		= "",
					xlab		= "Variance of TPM differential",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5)
									)   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
					)

		# Statistical test Kruskal Wallis)
		KW_result	<- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste(length(shared_genelist), " genes, P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")	
		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()
	}


	# Intraspecies diversity (with genes shared between stages)
	if(0){
		TPM_early_shared	<- na.omit(TPM_early)
		TPM_mid_shared		<- na.omit(TPM_mid)
		TPM_late1_shared	<- na.omit(TPM_late1)
		TPM_late2_shared	<- na.omit(TPM_late2)
		
		shared_genelist		<- intersect(rownames(TPM_early_shared), rownames(TPM_mid_shared))
		shared_genelist		<- intersect(shared_genelist, rownames(TPM_late1_shared))
		shared_genelist		<- intersect(shared_genelist, rownames(TPM_late2_shared))

		dist_early			<- Make_dist_array_wild(TPM_Kasasa_early[shared_genelist,], TPM_Oura_early[shared_genelist,])
		dist_mid			<- Make_dist_array_wild(TPM_Kasasa_mid[shared_genelist,], TPM_Oura_mid[shared_genelist,])
		dist_late1			<- Make_dist_array_wild(TPM_Kasasa_late1[shared_genelist,], TPM_Oura_late1[shared_genelist,])
		dist_late2			<- Make_dist_array_wild(TPM_Kasasa_late2[shared_genelist,], TPM_Oura_late2[shared_genelist,])

		outputpath	<- paste(outputdir_tmp2, "Box_wild_samegeneset.pdf", sep="")
		pdf(file = outputpath,	width = 7,		height = 7)
    	par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1, dist_late2,
									horizontal	= TRUE,
									ylim	= c(0, 0.15),
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					cex.axis	= 1,
									cex 	= 0.3,
									pch 	= 16,
									main 	= paste("After gene selection, Log, wild(only shared genes), genes with TPM=>",TPMCutOff, sep=""),
									las		= 1
								)

		par(new=T)

		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim	= c(0, 0.15),
					ylab		= "",
					xlab		= "Variance of TPM differential",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5))   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
				)

		axis(	side		= 4, 	
				at			= 1:4,  
    			labels		= c(paste( length(dist_early), " pairs\n(n = 4, 4)",sep=""),
    							paste( length(dist_mid), " pairs\n(n = 4, 4)",sep=""),
    							paste( length(dist_late1), " pairs\n(n = 4, 4)",sep=""),
    							paste( length(dist_late2), " pairs\n(n = 2, 2)",sep="")
    							),
    			tck			= 0,	
    			las			= 1,  	
   		 		mgp			= c(3.5,0.7,0), 
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)
		
		# Statistical test Kruskal Wallis)
		KW_result <- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste(length(shared_genelist), " genes, P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")

		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()

	}


	# Figures for explanation of calculation procedure
	if(1){
		outputpath	<- paste(outputdir_tmp2, "ExpressionLevels_st23.5_#2-1,2.pdf", sep="")
		pdf(file = outputpath, width = 7,	height = 7)
		par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 

		col_line	<- rep("dodgerblue", nrow(TPM_mid))
		for(k in 1:nrow(TPM_mid)){
			if(TPM_mid[k,"Ol_stage23.5_2.1"] < TPMCutOff || TPM_mid[k,"Ol_stage23.5_2.2"] < TPMCutOff){
				col_line[k]	<- "gray80"
			}
		}

		plot(	TPM_mid[,"Ol_stage23.5_2.1"], TPM_mid[,"Ol_stage23.5_2.2"],
				xlab	= "Expression level (st.23.5, twin #1-1)",
				ylab	= "Expression level (st.23.5, twin #1-2)",
				xlim	= c(0,4.5),
				ylim	= c(0,4.5),
				pch 	= 20,
				cex		= 0.2,
				col		= col_line,
				cex.lab = 1,
				main	= paste("After gene selection, Log", sep=""),
				las		= 1
			)
		lines(x=c(0,4.5), y=c(TPMCutOff,TPMCutOff), col="gray70")
		lines(x=c(TPMCutOff,TPMCutOff), y=c(0,4.5), col="gray70")

		dev.off()


		outputpath	<- paste(outputdir_tmp2, "ExpressionLevels_st15_#1-1,2.pdf", sep="")
		pdf(file = outputpath, width = 7,	height = 7)
		par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 

		col_line	<- rep("dodgerblue", nrow(TPM_early))
		for(k in 1:nrow(TPM_early)){
			if(TPM_early[k,"Ol_stage15_1.1"] < TPMCutOff || TPM_early[k,"Ol_stage15_1.2"] < TPMCutOff){
				col_line[k]	<- "gray80"
			}
		}

		plot(	TPM_early[,"Ol_stage15_1.1"], TPM_early[,"Ol_stage15_1.2"],
				xlab	= "Expression level (st.15, twin #1-1)",
				ylab	= "Expression level (st.15, twin #1-2)",
				xlim	= c(0,4.5),
				ylim	= c(0,4.5),
				pch 	= 20,
				cex		= 0.2,
				col		= col_line,
				cex.lab = 1,
				main	= paste("After gene selection, Log", sep=""),
				las		= 1
			)
		lines(x=c(0,4.5), y=c(TPMCutOff,TPMCutOff), col="gray70")
		lines(x=c(TPMCutOff,TPMCutOff), y=c(0,4.5), col="gray70")

		dev.off()



		get_expdifference	<- function(TPM_pair){
			TPM_pair	<- unlist(TPM_pair)
			if(TPM_pair[1] == 0 && TPM_pair[2] == 0){
				return(NA)
			}else{
				if(TPM_pair[1] < TPMCutOff || TPM_pair[2] < TPMCutOff){
					return(NA)
				}else{
					return((TPM_pair[1] - TPM_pair[2]))
				}
			}
		}

		tmp	<- apply(cbind(TPM_mid[,"Ol_stage23.5_2.1"],TPM_mid[,"Ol_stage23.5_2.2"]), 1, get_expdifference)

		outputpath	<- paste(outputdir_tmp2, "ExpressionDiff-and-Density_st23.5_#2-1,2.pdf", sep="")
		pdf(file = outputpath, width = 7,	height = 7)
		par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 

		plot(	density(tmp, na.rm=TRUE),
				xlim	= c(-1.8,1.8),
				ylim	= c(0,5),
				xlab	= "Expression Difference (st.23.5, betweem twin #1-1,2)",
				ylab	= "Density",
				main	= paste("After gene selection, Log, inbred, genes with TPM=>",TPMCutOff, sep=""),
				col		= "dodgerblue",
				lwd		= 2,
				las		= 1
			)
		#rug(tmp, col="deepskyblue")
		mtext(paste(length(na.omit(tmp)), " genes", sep=""), side = 1, line = 6, adj = 1)

		dev.off()
	}


	# Developmental stability (with threshold, by 1- Spearman's rho)
	if(1){	
		dist_early	<- Make_dist_array_spearman(TPM_early)
		dist_mid	<- Make_dist_array_spearman(TPM_mid)
		dist_late1	<- Make_dist_array_spearman(TPM_late1)
		dist_late2	<- Make_dist_array_spearman(TPM_late2)

		outputpath	<- paste(outputdir_tmp2, "Box_inbred_spearman_threshold.pdf", sep="")
		pdf(file = outputpath,	width = 7, height = 7)
   		par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1,dist_late2,
									horizontal	= TRUE,
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					#xaxt	= "n",
									ylim	= c(0, 0.1),
				 					cex.axis= 1,
									cex 	= 0.3,
									pch 	= 16,
									main	= paste("After gene selection, Log, inbred, genes with TPM=>",TPMCutOff, sep=""),
									las		= 1
								)

		axis(	side		= 4,
				at			= 1:4,
    			labels		= c(paste( length(dist_early), " pairs\n(n = ", ncol(TPM_early), ")",sep=""),
    							paste( length(dist_mid), " pairs\n(n = ", ncol(TPM_mid), ")",sep=""),
    							paste( length(dist_late1), " pairs\n(n = ", ncol(TPM_late1), ")",sep=""),
    							paste( length(dist_late2), " pairs\n(n = ", ncol(TPM_late2), ")",sep="")
    							),
    			tck			= 0,	
    			las			= 1,  		
   		 		mgp			= c(3.5,0.7,0), 
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)

		par(new=T)
			
		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))

		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xlim		= c(0, 0.1),
					xaxt		= "n",
					yaxt		= "n",
					ylab		= "",
					xlab		= "1- Spearman's rho",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5)
									)   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
					)

		# Statistical test Kruskal Wallis)
		KW_result	<- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")	
		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()
	}


	# Developmental stability (without threshold, by 1- Spearman's rho)
	if(1){	
		dist_early	<- Make_dist_array_spearman_0(TPM_early)
		dist_mid	<- Make_dist_array_spearman_0(TPM_mid)
		dist_late1	<- Make_dist_array_spearman_0(TPM_late1)
		dist_late2	<- Make_dist_array_spearman_0(TPM_late2)

		outputpath	<- paste(outputdir_tmp2, "Box_inbred_spearman_noThreshold.pdf", sep="")
		pdf(file = outputpath,	width = 7, height = 7)
   		par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1,dist_late2,
									horizontal	= TRUE,
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					#xaxt	= "n",
									ylim	= c(0, 0.1),
				 					cex.axis= 1,
									cex 	= 0.3,
									pch 	= 16,
									main	= paste("After gene selection, Log, inbred", sep=""),
									las		= 1
								)

		axis(	side		= 4,
				at			= 1:4,  
    			labels		= c(paste( length(dist_early), " pairs\n(n = ", ncol(TPM_early), ")",sep=""),
    							paste( length(dist_mid), " pairs\n(n = ", ncol(TPM_mid), ")",sep=""),
    							paste( length(dist_late1), " pairs\n(n = ", ncol(TPM_late1), ")",sep=""),
    							paste( length(dist_late2), " pairs\n(n = ", ncol(TPM_late2), ")",sep="")
    							),
    			tck			= 0,
    			las			= 1,  			）
   		 		mgp			= c(3.5,0.7,0), 		
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)

		par(new=T)
			
		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))

		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xlim		= c(0, 0.1),
					xaxt		= "n",
					yaxt		= "n",
					ylab		= "",
					xlab		= "1- Spearman's rho",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5)
									)   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
					)

		# Statistical test Kruskal Wallis)
		KW_result	<- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")	
		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()
	}


	# Intraspecies diversity (without threshold, by 1- Spearman's rho)
	if(1){
		dist_early			<- Make_dist_array_wild_spearman(TPM_Kasasa_early, TPM_Oura_early)
		dist_mid			<- Make_dist_array_wild_spearman(TPM_Kasasa_mid, TPM_Oura_mid)
		dist_late1			<- Make_dist_array_wild_spearman(TPM_Kasasa_late1, TPM_Oura_late1)
		dist_late2			<- Make_dist_array_wild_spearman(TPM_Kasasa_late2, TPM_Oura_late2)

		outputpath	<- paste(outputdir_tmp2, "Box_wild_spearman_NOthreshold.pdf", sep="")
		pdf(file = outputpath,	width = 7,		height = 7)
    	par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1, dist_late2,
									horizontal	= TRUE,
									ylim	= c(0, 0.15),
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					cex.axis	= 1,
									cex 	= 0.3,
									pch 	= 16,
									main 	= paste("After gene selection, Log, wild, genes with TPM=>",TPMCutOff, sep=""),
									las		= 1
								)

		par(new=T)

		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim	= c(0, 0.15),
					ylab		= "",
					xlab		= "1 - Spearman's rho",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5))   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
				)

		axis(	side		= 4, 	
				at			= 1:4, 
    			labels		= c(paste( length(dist_early), " pairs\n(n = 4, 4)",sep=""),
    							paste( length(dist_mid), " pairs\n(n = 4, 4)",sep=""),
    							paste( length(dist_late1), " pairs\n(n = 4, 4)",sep=""),
    							paste( length(dist_late2), " pairs\n(n = 2, 2)",sep="")
    							),
    			tck			= 0,		
    			las			= 1,  			
   		 		mgp			= c(3.5,0.7,0), 	
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)
		
		# Statistical test Kruskal Wallis)
		KW_result <- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")

		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()

	}




	rm(box_position)
	rm(dist)
	rm(stage)
	rm(DistStage)
	rm(dist_early)
	rm(dist_mid)
	rm(dist_late1)
	rm(dist_late2)
	rm(dist_early_tech)
	rm(dist_mid_tech)
	rm(dist_late1_tech)
	rm(dist_late2_tech)

	rm(TPM_early_F)
	rm(TPM_mid_F)
	rm(TPM_late1_F)
	rm(TPM_late2_F)
	rm(TPM_early_M)
	rm(TPM_mid_M)
	rm(TPM_late1_M)
	rm(TPM_late2_M)

	rm(TPM_early_shared)
	rm(TPM_mid_shared)
	rm(TPM_late1_shared)
	rm(TPM_late2_shared)
	rm(shared_genelist)

	rm(tmp)

	gc()
	gc()
}


##########################################################################################################################################
#### Calculate each of developmental stability, intraspecies diversity, and tech error; then correct for expression level dependence
##########################################################################################################################################

# Calculation
if(1){
	stages	<- c("early", "mid", "late1", "late2")
	for(i in 1:4){
		stage			<- stages[i]

		tmp_inbred		<- na.omit(get(paste("TPM_", stage, sep="")))	
		tmp_tech		<- get(paste("TPM_", stage, "_pool", sep=""))[rownames(tmp_inbred),]
		tmp_Kasasa		<- get(paste("TPM_Kasasa_", stage, sep=""))[rownames(tmp_inbred),]
		tmp_Oura		<- get(paste("TPM_Oura_", stage, sep=""))[rownames(tmp_inbred),]



		## Sort the genes by expression level and obtain the running median of developmental stability. Then correct for gene expression level dependence
			tmp_flu_nocorrection	<- apply(tmp_inbred, 1, get_fluctuation_log)
			tmp_flu_nocorrection	<- na.omit(tmp_flu_nocorrection)								
			tmp_explevel			<- apply(tmp_inbred[names(tmp_flu_nocorrection),], 1, mean)	
			tmp_explevel			<- tmp_explevel[order(tmp_explevel, decreasing=TRUE)]			
			tmp_flu_nocorrection	<- tmp_flu_nocorrection[names(tmp_explevel)]				

			RM_flu			<- numeric(length(tmp_flu_nocorrection))
			names(RM_flu)	<- names(tmp_flu_nocorrection)
			for(k in 1:length(tmp_flu_nocorrection)){
				
				if(k == 1 || k == length(tmp_flu_nocorrection)){
					RM_flu[k]	<- tmp_flu_nocorrection[k]
				}else{
				
					if(k <= 250){
						RM_flu[k]	<- median(tmp_flu_nocorrection[1:(2*k-1)], na.rm=TRUE)
					}else{
						if(k <= (length(tmp_flu_nocorrection) - 250)){
							RM_flu[k]	<- median(tmp_flu_nocorrection[(k-250):(k+250)], na.rm=TRUE)
						}else{
							RM_flu[k]	<- median(tmp_flu_nocorrection[(2*k-length(tmp_flu_nocorrection)):length(tmp_flu_nocorrection)], na.rm=TRUE)
						}						
					}					
				}	
			}

			tmp_flu_corrected		<- tmp_flu_nocorrection - RM_flu



		## Sort the genes by expression level and obtain the running median of intraspecies diversity. Then correct for gene expression level dependence
			tmp_microevo_nocorrection_new	<- apply(cbind(tmp_Kasasa,tmp_Oura), 1, get_microevo_log2)
			tmp_microevo_nocorrection_new	<- tmp_microevo_nocorrection_new[names(tmp_flu_nocorrection)]
			tmp_microevo_nocorrection_new	<- tmp_microevo_nocorrection_new[names(tmp_explevel)]	

			RM_evo_new			<- numeric(length(tmp_microevo_nocorrection_new))
			names(RM_evo_new)	<- names(tmp_microevo_nocorrection_new)
			for(k in 1:length(tmp_microevo_nocorrection_new)){
				
				if(k == 1 || k == length(tmp_microevo_nocorrection_new)){
					RM_evo_new[k]	<- tmp_microevo_nocorrection_new[k]
				}else{
				
					if(k <= 250){
						RM_evo_new[k]	<- median(tmp_microevo_nocorrection_new[1:(2*k-1)], na.rm=TRUE)
					}else{
						if(k <= (length(tmp_microevo_nocorrection_new) - 250)){
							RM_evo_new[k]	<- median(tmp_microevo_nocorrection_new[(k-250):(k+250)], na.rm=TRUE)
						}else{
							RM_evo_new[k]	<- median(tmp_microevo_nocorrection_new[(2*k-length(tmp_microevo_nocorrection_new)):length(tmp_microevo_nocorrection_new)], na.rm=TRUE)
						}						
					}					
				}	
			}

			tmp_microevo_nocorrection_new			<- na.omit(tmp_microevo_nocorrection_new)
			RM_evo_new								<- RM_evo_new[names(tmp_microevo_nocorrection_new)]
			tmp_microevo_corrected_new				<- tmp_microevo_nocorrection_new - RM_evo_new


		## Sort the genes by expression level and obtain the running median of technical error. Then correct for gene expression level dependence
			tmp_techerror_nocorrection	<- apply(tmp_tech, 1, get_techerror_log)
			tmp_techerror_nocorrection	<- tmp_techerror_nocorrection[names(tmp_flu_nocorrection)]	
			tmp_techerror_nocorrection	<- tmp_techerror_nocorrection[names(tmp_explevel)]		

			RM_tech			<- numeric(length(tmp_techerror_nocorrection))
			names(RM_tech)	<- names(tmp_techerror_nocorrection)
			for(k in 1:length(tmp_techerror_nocorrection)){
				
				if(k == 1 || k == length(tmp_techerror_nocorrection)){
					RM_tech[k]	<- tmp_techerror_nocorrection[k]
				}else{
				
					if(k <= 250){
						RM_tech[k]	<- median(tmp_techerror_nocorrection[1:(2*k-1)], na.rm=TRUE)
					}else{
						if(k <= (length(tmp_techerror_nocorrection) - 250)){
							RM_tech[k]	<- median(tmp_techerror_nocorrection[(k-250):(k+250)], na.rm=TRUE)
						}else{
							RM_tech[k]	<- median(tmp_techerror_nocorrection[(2*k-length(tmp_techerror_nocorrection)):length(tmp_techerror_nocorrection)], na.rm=TRUE)
						}						
					}					
				}	
			}																		

			tmp_techerror_nocorrection	<- na.omit(tmp_techerror_nocorrection)	
			RM_tech						<- RM_tech[names(tmp_techerror_nocorrection)]	
			tmp_techerror_corrected		<- tmp_techerror_nocorrection - RM_tech


	assign(paste("flu_nocorrection_", stage, sep=""), tmp_flu_nocorrection)
	assign(paste("RM_flu_", stage, sep=""), RM_flu)
	assign(paste("flu_RMcorrection_", stage, sep=""), tmp_flu_corrected)

	assign(paste("microevo_nocorrection_conventional_", stage, sep=""), tmp_microevo_nocorrection_conventional)
	assign(paste("RM_evo_conventional_", stage, sep=""), RM_evo_conventional)
	assign(paste("microevo_RMcorrection_conventional_", stage, sep=""), tmp_microevo_corrected_conventional)

	assign(paste("microevo_nocorrection_new_", stage, sep=""), tmp_microevo_nocorrection_new)
	assign(paste("RM_evo_new_", stage, sep=""), RM_evo_new)
	assign(paste("microevo_RMcorrection_new_", stage, sep=""), tmp_microevo_corrected_new)

	assign(paste("techerror_nocorrection_", stage, sep=""), tmp_techerror_nocorrection)
	assign(paste("RM_tech_", stage, sep=""), RM_tech)
	assign(paste("techerror_RMcorrection_", stage, sep=""), tmp_techerror_corrected)

	}
}

## Visualization
if(1){
	outputdir_tmp2	<- paste(outputdir_tmp, "fluctuation,microevo,explevel,techerror_AfterGeneSelection/", sep="")
	dir.create(outputdir_tmp2)

	stages	<- c("early", "mid", "late1", "late2")
	for(i in 1:4){
		stage										<- stages[i]
		tmp_inbred									<- na.omit(get(paste("TPM_", stage, sep="")))

		## 揺らぎ関係
		tmp_flu_nocorrection					<- get(paste("flu_nocorrection_", stage, sep=""))
		RM_flu									<- get(paste("RM_flu_", stage, sep=""))
		tmp_flu_corrected						<- get(paste("flu_RMcorrection_", stage, sep=""))
		outputpath								<- paste(outputdir_tmp2, "RM_fluctuation_", stage, ".txt", sep="")
		write.table(data.frame(cbind(gene_names[names(tmp_flu_corrected)], tmp_flu_corrected)), file = outputpath, quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE )


		## 小進化変化関係
		tmp_microevo_nocorrection_new			<- get(paste("microevo_nocorrection_new_", stage, sep=""))
		RM_evo_new								<- get(paste("RM_evo_new_", stage, sep=""))
		tmp_microevo_corrected_new				<- get(paste("microevo_RMcorrection_new_", stage, sep=""))
		outputpath								<- paste(outputdir_tmp2, "RM_microevo_", stage, ".txt", sep="")
		write.table(data.frame(cbind(gene_names[names(tmp_microevo_corrected_new)], tmp_microevo_corrected_new)), file = outputpath, quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE )


		## tech error関係
		tmp_techerror_nocorrection				<- get(paste("techerror_nocorrection_", stage, sep=""))	
		RM_tech									<- get(paste("RM_tech_", stage, sep=""))
		tmp_techerror_corrected					<- get(paste("techerror_RMcorrection_", stage, sep=""))
		
		## 発現量を求める。
		tmp_explevel							<- apply(tmp_inbred[names(tmp_flu_nocorrection),], 1, mean)	


		## 描画　修正前の揺らぎと発現量、RMと発現量を重ね書き
		if(1){
			outputpath	<- paste(outputdir_tmp2, "NoRMCorrection_fluctuation-explevel_RM-explevel_", stage, ".pdf", sep="")
			pdf(file = outputpath, width = 7,	height = 7)
			par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 

			plot(	tmp_explevel[names(tmp_flu_nocorrection)], tmp_flu_nocorrection,  
					xlim	= c(0, 5), 
					ylim	= c(0, 1.7), 
					#xaxt	= "n",
					#yaxt	= "n",
					#log		= "xy",
					xlab	= "Expression level",
					ylab	= "Fluctuation (no RM correction), RM",
					pch 	= 20,
					cex		= 0.2,
					col		= "dodgerblue",
					cex.lab = 1,
					main	= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
					las		=1
				)

			par(new = TRUE)
			
			plot(	tmp_explevel[names(tmp_flu_nocorrection)], RM_flu, 
					xlim	= c(0, 5), 
					ylim	= c(0, 1.7), 
					xaxt	= "n",
					yaxt	= "n",
					#log		= "xy",
					xlab	= "",
					ylab	= "",
					pch 	= 20,
					cex		= 0.2,
					col		= "orangered",
					cex.lab = 0.8
				)

			mtext(paste(length(na.omit(tmp_flu_nocorrection)), " genes", sep=""), side = 1, line = 6, adj = 1)

			cortest_list	<- cor.test(tmp_explevel[names(tmp_flu_nocorrection)], tmp_flu_nocorrection, method="spearman")
			mtext(paste("fluctuation(no RM correction), explevel rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
			
			dev.off()
		}


		## 描画　修正前の小進化応答と発現量、RMと発現量を重ね書き
		if(1){
			outputpath	<- paste(outputdir_tmp2, "NoRMCorrection_microevo(new)-explevel_RM-explevel_", stage, ".pdf", sep="")
			pdf(file = outputpath, width = 7,	height = 7)
			par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 

			plot(	tmp_explevel[names(tmp_microevo_nocorrection_new)], tmp_microevo_nocorrection_new, 
					xlim	= c(0, 5), 
					ylim	= c(0, 2), 
					#xaxt	= "n",
					#yaxt	= "n",
					#log		= "xy",
					xlab	= "Expression level",
					ylab	= "Microevolutionary diversity (no RM correction, |i-j|->mean), RM",
					pch 	= 20,
					cex		= 0.2,
					col		= "dodgerblue",
					cex.lab = 0.8,
					main	= paste("After gene selection, Log, meanTPM->",TPMCutOff, ", ", stage, sep=""),
					las		=1
				)

			par(new = TRUE)
			
			plot(	tmp_explevel[names(tmp_microevo_nocorrection_new)], RM_evo_new, 
					xlim	= c(0, 5), 
					ylim	= c(0, 2), 
					xaxt	= "n",
					yaxt	= "n",
					#log		= "xy",
					xlab	= "",
					ylab	= "",
					pch 	= 20,
					cex		= 0.2,
					col		= "orangered",
					cex.lab = 0.8
				)

			mtext(paste(length(na.omit(tmp_microevo_nocorrection_new)), " genes", sep=""), side = 1, line = 6, adj = 1)

			cortest_list	<- cor.test(tmp_explevel[names(tmp_microevo_nocorrection_new)], tmp_microevo_nocorrection_new, method="spearman")
			mtext(paste("microevo(no RM correction, new)-explevel rho = ", signif(cortest_list$estimate, 2), " (P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
			
			dev.off()
		}


		## 描画　修正前のtech errorと発現量、RMと発現量を重ね書き
		if(1){
			outputpath	<- paste(outputdir_tmp2, "NoRMCorrection_techerror-explevel_RM-explevel_", stage, ".pdf", sep="")
			pdf(file = outputpath, width = 7,	height = 7)
			par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 

			plot(	tmp_explevel[names(tmp_techerror_nocorrection)], tmp_techerror_nocorrection, 
					xlim	= c(0, 5), 
					ylim	= c(0, 1.6), 
					#xaxt	= "n",
					#yaxt	= "n",
					#log		= "xy",
					xlab	= "Expression level",
					ylab	= "Tech error (no RM correction), RM",
					pch 	= 20,
					cex		= 0.2,
					col		= "dodgerblue",
					cex.lab = 1,
					main	= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
					las		=1
				)

			par(new = TRUE)
			
			plot(	tmp_explevel[names(tmp_techerror_nocorrection)], RM_tech, 
					xlim	= c(0, 5), 
					ylim	= c(0, 1.6), 
					xaxt	= "n",
					yaxt	= "n",
					#log		= "xy",
					xlab	= "",
					ylab	= "",
					pch 	= 20,
					cex		= 0.2,
					col		= "orangered",
					cex.lab = 0.8
				)

			mtext(paste(length(na.omit(tmp_techerror_nocorrection)), " genes", sep=""), side = 1, line = 6, adj = 1)

			cortest_list	<- cor.test(tmp_explevel[names(tmp_techerror_nocorrection)], tmp_techerror_nocorrection, method="spearman")
			mtext(paste("techerror(no RM correction)-explevel rho = ", signif(cortest_list$estimate, 2), " (P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
			
			dev.off()
		}

		## 描画　修正前の揺らぎと修正前の小進化変化、修正前の揺らぎと修正後のtech error
		if(1){
			tmp_gene	<- intersect(names(tmp_flu_nocorrection), names(tmp_microevo_nocorrection_new))
			tmp_gene	<- intersect(tmp_gene, names(tmp_techerror_nocorrection))

			outputpath	<- paste(outputdir_tmp2, "NORMCorrection_fluctuation-microevo(new)_fluctuation-techerror_", stage, ".pdf", sep="")
			pdf(file = outputpath, width = 7,	height = 7)
			par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 
																												
			#color  = densCols(tmp_flu_corrected[tmp_gene], tmp_microevo_corrected_new[tmp_gene], colramp = colorRampPalette(c("#01579B","#0277BD","#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
			color  = densCols(tmp_flu_nocorrection[tmp_gene], tmp_microevo_nocorrection_new[tmp_gene], colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
			plot(	tmp_flu_nocorrection[tmp_gene], tmp_microevo_nocorrection_new[tmp_gene], 
					xlim	= c(-0.2, 0.8), 
					ylim	= c(-0.2, 1.2), 
					#xaxt	= "n",
					#yaxt	= "n",
					#log		= "xy",
					xlab	= "Fluctuation (NO correction)",
					ylab	= "Microevolutionary diversity (No correction, new), Tech error (No correction)",
					pch 	= 20,
					cex		= 0.2,	#0.08,
					col		= color, #rgb(0/255, 191/255, 255/255, alpha=0.1),	#"dodgerblue"
					cex.lab = 1,
					main	= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
					las		=1
				)
			
			par(new=T)


			#color  = densCols(tmp_flu_corrected[tmp_gene], tmp_techerror_corrected[tmp_gene], colramp = colorRampPalette(c("khaki3","khaki1")), nbin = 3000)
			color  = densCols(tmp_flu_nocorrection[tmp_gene], tmp_techerror_nocorrection[tmp_gene], colramp = colorRampPalette(c("gold3","gold1")), nbin = 3000)
			plot(	tmp_flu_nocorrection[tmp_gene], tmp_techerror_nocorrection[tmp_gene], 
					xlim	= c(-0.2, 0.8), 
					ylim	= c(-0.2, 1.2), 
					#xaxt	= "n",
					#yaxt	= "n",
					#log		= "xy",
					xlab	= "",
					ylab	= "",
					pch 	= 20,
					cex		= 0.05,	#0.08,
					col		= color, #rgb(204/255, 204/255, 204/255, alpha=0.5),
					cex.lab = 1,
					las		=1
				)
		
			mtext(paste(nrow(na.omit(cbind(tmp_flu_nocorrection[tmp_gene], tmp_microevo_nocorrection_new[tmp_gene]))), " genes", sep=""), side = 1, line = 6, adj = 1)

			cortest_list	<- cor.test(tmp_flu_nocorrection[tmp_gene], tmp_microevo_nocorrection_new[tmp_gene], method="spearman")
			cortest_list2	<- cor.test(tmp_flu_nocorrection[tmp_gene], tmp_techerror_nocorrection[tmp_gene], method="spearman")
			mtext(paste("flu-evo rho = ", signif(cortest_list$estimate, 2), " (P = ", signif(cortest_list$p.value, 2),"); flu-tech rho = ", signif(cortest_list2$estimate, 2), " (P = ", signif(cortest_list2$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
			
			dev.off()


		}

		## 描画　修正後の揺らぎと修正後の小進化変化、修正後の揺らぎと修正後のtech error
		if(1){
			tmp_gene	<- intersect(names(tmp_flu_corrected), names(tmp_microevo_corrected_new))
			tmp_gene	<- intersect(tmp_gene, names(tmp_techerror_corrected))

			outputpath	<- paste(outputdir_tmp2, "RMCorrection_fluctuation-microevo(new)_fluctuation-techerror_", stage, ".pdf", sep="")
			pdf(file = outputpath, width = 7,	height = 7)
			par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 
																												
			#color  = densCols(tmp_flu_corrected[tmp_gene], tmp_microevo_corrected_new[tmp_gene], colramp = colorRampPalette(c("#01579B","#0277BD","#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
			color  = densCols(tmp_flu_corrected[tmp_gene], tmp_microevo_corrected_new[tmp_gene], colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
			plot(	tmp_flu_corrected[tmp_gene], tmp_microevo_corrected_new[tmp_gene], 
					xlim	= c(-0.2, 0.8), 
					ylim	= c(-0.2, 1.2), 
					#xaxt	= "n",
					#yaxt	= "n",
					#log		= "xy",
					xlab	= "Fluctuation (RM correction)",
					ylab	= "Microevolutionary diversity (RM correction, new), Tech error (RM correction)",
					pch 	= 20,
					cex		= 0.2,	#0.08,
					col		= color, #rgb(0/255, 191/255, 255/255, alpha=0.1),	#"dodgerblue"
					cex.lab = 1,
					main	= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
					las		=1
				)
			
			par(new=T)


			#color  = densCols(tmp_flu_corrected[tmp_gene], tmp_techerror_corrected[tmp_gene], colramp = colorRampPalette(c("khaki3","khaki1")), nbin = 3000)
			color  = densCols(tmp_flu_corrected[tmp_gene], tmp_techerror_corrected[tmp_gene], colramp = colorRampPalette(c("gold3","gold1")), nbin = 3000)
			plot(	tmp_flu_corrected[tmp_gene], tmp_techerror_corrected[tmp_gene], 
					xlim	= c(-0.2, 0.8), 
					ylim	= c(-0.2, 1.2), 
					#xaxt	= "n",
					#yaxt	= "n",
					#log		= "xy",
					xlab	= "",
					ylab	= "",
					pch 	= 20,
					cex		= 0.05,	#0.08,
					col		= color, #rgb(204/255, 204/255, 204/255, alpha=0.5),
					cex.lab = 1,
					las		=1
				)
		
			mtext(paste(nrow(na.omit(cbind(tmp_flu_corrected[tmp_gene], tmp_microevo_corrected_new[tmp_gene]))), " genes", sep=""), side = 1, line = 6, adj = 1)

			cortest_list	<- cor.test(tmp_flu_corrected[tmp_gene], tmp_microevo_corrected_new[tmp_gene], method="spearman")
			cortest_list2	<- cor.test(tmp_flu_corrected[tmp_gene], tmp_techerror_corrected[tmp_gene], method="spearman")
			mtext(paste("flu-evo rho = ", signif(cortest_list$estimate, 2), " (P = ", signif(cortest_list$p.value, 2),"); flu-tech rho = ", signif(cortest_list2$estimate, 2), " (P = ", signif(cortest_list2$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
			
			dev.off()


		}


		## 描画　修正後の揺らぎと発現量、修正後のtech errorと発現量
		if(1){
			tmp_gene	<- intersect(names(tmp_flu_corrected), names(tmp_techerror_corrected))

			outputpath	<- paste(outputdir_tmp2, "RMCorrection_explevel-fluctuation_explevel-techerror_", stage, ".pdf", sep="")
			pdf(file = outputpath, width = 7,	height = 7)
			par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 

			#color  = densCols(tmp_explevel[names(tmp_flu_corrected)], tmp_flu_corrected, colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
			plot(	tmp_explevel[names(tmp_flu_corrected)], tmp_flu_corrected, 
					xlim	= c(0, 5), 
					ylim	= c(-0.2, 0.8), 
					#xaxt	= "n",
					#yaxt	= "n",
					#log		= "xy",
					xlab	= "Expression level",
					ylab	= "Fluctuation (RM correction), Tech error (RM correction)",
					pch 	= 20,
					cex		= 0.1, #0.08,
					col		= "dodgerblue",
					cex.lab = 1,
					main	= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
					las		=1
				)
			if(0){
				par(new=T)

				plot(	tmp_explevel[tmp_gene], tmp_techerror_corrected[tmp_gene], 
					xlim	= c(0, 5), 
					ylim	= c(-0.2, 0.8), 
					xaxt	= "n",
					yaxt	= "n",
					#log		= "xy",
					xlab	= "",
					ylab	= "",
					pch 	= 20,
					cex		= 0.1, #0.08,
					col		= "gray80", #rgb(127/255, 127/255, 127/255, alpha=0.3),
					cex.lab = 1,
					las		=1
					)
			}

			mtext(paste(nrow(na.omit(cbind(tmp_explevel[tmp_gene], tmp_flu_corrected[tmp_gene]))), " genes", sep=""), side = 1, line = 6, adj = 1)

			cortest_list	<- cor.test(tmp_explevel[tmp_gene], tmp_flu_corrected[tmp_gene], method="spearman")
			mtext(paste("fluctuation(RM correction only)-exp level rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
			
			dev.off()
		}

		## 描画　修正後の小進化応答と発現量
		if(1){
			outputpath	<- paste(outputdir_tmp2, "RMCorrection_explevel-microevo(new)_", stage, ".pdf", sep="")
			pdf(file = outputpath, width = 7,	height = 7)
			par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 

			#color  = densCols(tmp_explevel[names(tmp_microevo_corrected_new)], tmp_microevo_corrected_new, colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
			plot(	tmp_explevel[names(tmp_microevo_corrected_new)], tmp_microevo_corrected_new, 
					xlim	= c(0, 5), 
					ylim	= c(-0.2, 1.2), 
					#xaxt	= "n",
					#yaxt	= "n",
					#log		= "xy",
					xlab	= "Expression level",
					ylab	= "Microevolutionary diversity (RM correction, new)",
					pch 	= 20,
					cex		= 0.1, #0.08,
					col		= "dodgerblue",
					cex.lab = 1,
					main	= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
					las		=1
				)
			
			mtext(paste(nrow(na.omit(cbind(tmp_explevel[names(tmp_microevo_corrected_new)], tmp_microevo_corrected_new))), " genes", sep=""), side = 1, line = 6, adj = 1)

			cortest_list	<- cor.test(tmp_explevel[names(tmp_microevo_corrected_new)], tmp_microevo_corrected_new, method="spearman")
			mtext(paste("microevo(RM correction)-exp level rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
			
			dev.off()
		}

	}
}		


# Calculation of intraspecies diversity between inbred-Kasasa、inbred-Oura
if(0){
	stages	<- c("early", "mid", "late1", "late2")
	for(i in 1:4){
		stage			<- stages[i]

		tmp_inbred		<- na.omit(get(paste("TPM_", stage, sep="")))	
		tmp_Kasasa		<- get(paste("TPM_Kasasa_", stage, sep=""))[rownames(tmp_inbred),]
		tmp_Oura		<- get(paste("TPM_Oura_", stage, sep=""))[rownames(tmp_inbred),]


        ## inbred-Kasasa
            tmp_microevo_nocorrection_Kasasa    <- apply(cbind(tmp_Kasasa,tmp_inbred), 1, get_microevo_log2)
            tmp_microevo_nocorrection_Kasasa    <- tmp_microevo_nocorrection_Kasasa[names(tmp_flu_nocorrection)]  
            tmp_microevo_nocorrection_Kasasa    <- tmp_microevo_nocorrection_Kasasa[names(tmp_explevel)]     

            RM_evo_Kasasa           <- numeric(length(tmp_microevo_nocorrection_Kasasa))
            names(RM_evo_Kasasa)    <- names(tmp_microevo_nocorrection_Kasasa)
            for(k in 1:length(tmp_microevo_nocorrection_Kasasa)){
        
                if(k == 1 || k == length(tmp_microevo_nocorrection_Kasasa)){
                    RM_evo_Kasasa[k]    <- tmp_microevo_nocorrection_Kasasa[k]
                }else{
               
                    if(k <= 250){
                        RM_evo_Kasasa[k]   <- median(tmp_microevo_nocorrection_Kasasa[1:(2*k-1)], na.rm=TRUE)
                    }else{
                        if(k <= (length(tmp_microevo_nocorrection_Kasasa) - 250)){
                            RM_evo_Kasasa[k]   <- median(tmp_microevo_nocorrection_Kasasa[(k-250):(k+250)], na.rm=TRUE)
                        }else{
                            RM_evo_Kasasa[k]   <- median(tmp_microevo_nocorrection_Kasasa[(2*k-length(tmp_microevo_nocorrection_Kasasa)):length(tmp_microevo_nocorrection_Kasasa)], na.rm=TRUE)
                        }                       
                    }                   
                }   
            }

            tmp_microevo_nocorrection_Kasasa           <- na.omit(tmp_microevo_nocorrection_Kasasa)
            RM_evo_Kasasa                              <- RM_evo_Kasasa[names(tmp_microevo_nocorrection_Kasasa)]
            tmp_microevo_corrected_Kasasa              <- tmp_microevo_nocorrection_Kasasa - RM_evo_Kasasa

    assign(paste("microevo_nocorrection_Kasasa_", stage, sep=""), tmp_microevo_nocorrection_Kasasa)
    assign(paste("RM_evo_Kasasa_", stage, sep=""), RM_evo_Kasasa)
    assign(paste("microevo_RMcorrection_Kasasa_", stage, sep=""), tmp_microevo_corrected_Kasasa)

        ## inbred-Kasasa
            tmp_microevo_nocorrection_Oura    <- apply(cbind(tmp_Oura,tmp_inbred), 1, get_microevo_log2)
            tmp_microevo_nocorrection_Oura    <- tmp_microevo_nocorrection_Oura[names(tmp_flu_nocorrection)]    
            tmp_microevo_nocorrection_Oura    <- tmp_microevo_nocorrection_Oura[names(tmp_explevel)]     

            RM_evo_Oura           <- numeric(length(tmp_microevo_nocorrection_Oura))
            names(RM_evo_Oura)    <- names(tmp_microevo_nocorrection_Oura)
            for(k in 1:length(tmp_microevo_nocorrection_Oura)){
               
                if(k == 1 || k == length(tmp_microevo_nocorrection_Oura)){
                    RM_evo_Oura[k]    <- tmp_microevo_nocorrection_Oura[k]
                }else{
                
                    if(k <= 250){
                        RM_evo_Oura[k]   <- median(tmp_microevo_nocorrection_Oura[1:(2*k-1)], na.rm=TRUE)
                    }else{
                        if(k <= (length(tmp_microevo_nocorrection_Oura) - 250)){
                            RM_evo_Oura[k]   <- median(tmp_microevo_nocorrection_Oura[(k-250):(k+250)], na.rm=TRUE)
                        }else{
                            RM_evo_Oura[k]   <- median(tmp_microevo_nocorrection_Oura[(2*k-length(tmp_microevo_nocorrection_Oura)):length(tmp_microevo_nocorrection_Oura)], na.rm=TRUE)
                        }                       
                    }                   
                }   
            }

            tmp_microevo_nocorrection_Oura           <- na.omit(tmp_microevo_nocorrection_Oura)
            RM_evo_Oura                              <- RM_evo_Oura[names(tmp_microevo_nocorrection_Oura)]
            tmp_microevo_corrected_Oura              <- tmp_microevo_nocorrection_Oura - RM_evo_Oura

    assign(paste("microevo_nocorrection_Oura_", stage, sep=""), tmp_microevo_nocorrection_Oura)
    assign(paste("RM_evo_Oura_", stage, sep=""), RM_evo_Oura)
    assign(paste("microevo_RMcorrection_Oura_", stage, sep=""), tmp_microevo_corrected_Oura)
	}
}

## Visualization
if(0){
	outputdir_tmp2	<- paste(outputdir_tmp, "fluctuation,microevo,explevel,techerror_AfterGeneSelection/", sep="")
	dir.create(outputdir_tmp2)

	stages	<- c("early", "mid", "late1", "late2")
	for(i in 1:4){
		stage										<- stages[i]
		tmp_inbred									<- na.omit(get(paste("TPM_", stage, sep="")))

		tmp_flu_nocorrection					<- get(paste("flu_nocorrection_", stage, sep=""))
		RM_flu									<- get(paste("RM_flu_", stage, sep=""))
		tmp_flu_corrected						<- get(paste("flu_RMcorrection_", stage, sep=""))

        tmp_microevo_nocorrection_Kasasa           <- get(paste("microevo_nocorrection_Kasasa_", stage, sep=""))
        RM_evo_Kasasa                              <- get(paste("RM_evo_Kasasa_", stage, sep=""))
        tmp_microevo_corrected_Kasasa              <- get(paste("microevo_RMcorrection_Kasasa_", stage, sep=""))

        tmp_microevo_nocorrection_Oura           <- get(paste("microevo_nocorrection_Oura_", stage, sep=""))
        RM_evo_Oura                              <- get(paste("RM_evo_Oura_", stage, sep=""))
        tmp_microevo_corrected_Oura              <- get(paste("microevo_RMcorrection_Oura_", stage, sep=""))

		tmp_techerror_nocorrection				<- get(paste("techerror_nocorrection_", stage, sep=""))	
		RM_tech									<- get(paste("RM_tech_", stage, sep=""))
		tmp_techerror_corrected					<- get(paste("techerror_RMcorrection_", stage, sep=""))

		tmp_explevel							<- apply(tmp_inbred[names(tmp_flu_nocorrection),], 1, mean)	


        if(1){
            outputpath  <- paste(outputdir_tmp2, "NoRMCorrection_microevo(inbred-Kasasa)-explevel_RM-explevel_", stage, ".pdf", sep="")
            pdf(file = outputpath, width = 7,   height = 7)
            par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 

            plot(   tmp_explevel[names(tmp_microevo_nocorrection_Kasasa)], tmp_microevo_nocorrection_Kasasa, 
                    xlim    = c(0, 5), 
                    ylim    = c(0, 2), 
                    #xaxt   = "n",
                    #yaxt   = "n",
                    #log        = "xy",
                    xlab    = "Expression level",
                    ylab    = "Microevolutionary diversity (no RM correction, inbred-Kasasa), RM",
                    pch     = 20,
                    cex     = 0.2,
                    col     = "dodgerblue",
                    cex.lab = 0.8,
                    main    = paste("After gene selection, Log, meanTPM->",TPMCutOff, ", ", stage, sep=""),
                    las     =1
                )

            par(new = TRUE)
            
            plot(   tmp_explevel[names(tmp_microevo_nocorrection_Kasasa)], RM_evo_Kasasa, 
                    xlim    = c(0, 5), 
                    ylim    = c(0, 2), 
                    xaxt    = "n",
                    yaxt    = "n",
                    #log        = "xy",
                    xlab    = "",
                    ylab    = "",
                    pch     = 20,
                    cex     = 0.2,
                    col     = "orangered",
                    cex.lab = 0.8
                )

            mtext(paste(length(na.omit(tmp_microevo_nocorrection_Kasasa)), " genes", sep=""), side = 1, line = 6, adj = 1)

            cortest_list    <- cor.test(tmp_explevel[names(tmp_microevo_nocorrection_Kasasa)], tmp_microevo_nocorrection_Kasasa, method="spearman")
            mtext(paste("microevo(no RM correction, inbred-Kasasa)-explevel rho = ", signif(cortest_list$estimate, 2), " (P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
            
            dev.off()
        }

       	if(1){
            outputpath  <- paste(outputdir_tmp2, "NoRMCorrection_microevo(inbred-Oura)-explevel_RM-explevel_", stage, ".pdf", sep="")
            pdf(file = outputpath, width = 7,   height = 7)
            par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 

            plot(   tmp_explevel[names(tmp_microevo_nocorrection_Oura)], tmp_microevo_nocorrection_Oura, 
                    xlim    = c(0, 5), 
                    ylim    = c(0, 2), 
                    #xaxt   = "n",
                    #yaxt   = "n",
                    #log        = "xy",
                    xlab    = "Expression level",
                    ylab    = "Microevolutionary diversity (no RM correction, inbred-Oura), RM",
                    pch     = 20,
                    cex     = 0.2,
                    col     = "dodgerblue",
                    cex.lab = 0.8,
                    main    = paste("After gene selection, Log, meanTPM->",TPMCutOff, ", ", stage, sep=""),
                    las     =1
                )

            par(new = TRUE)
            
            plot(   tmp_explevel[names(tmp_microevo_nocorrection_Oura)], RM_evo_Oura, 
                    xlim    = c(0, 5), 
                    ylim    = c(0, 2), 
                    xaxt    = "n",
                    yaxt    = "n",
                    #log        = "xy",
                    xlab    = "",
                    ylab    = "",
                    pch     = 20,
                    cex     = 0.2,
                    col     = "orangered",
                    cex.lab = 0.8
                )

            mtext(paste(length(na.omit(tmp_microevo_nocorrection_Oura)), " genes", sep=""), side = 1, line = 6, adj = 1)

            cortest_list    <- cor.test(tmp_explevel[names(tmp_microevo_nocorrection_Oura)], tmp_microevo_nocorrection_Oura, method="spearman")
            mtext(paste("microevo(no RM correction, inbred-Oura)-explevel rho = ", signif(cortest_list$estimate, 2), " (P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
            
            dev.off()
        }


         if(1){
            tmp_gene    <- intersect(names(tmp_flu_corrected), names(tmp_microevo_corrected_Kasasa))

            outputpath  <- paste(outputdir_tmp2, "RMCorrection_fluctuation-microevo(inbred-Kasasa)_", stage, ".pdf", sep="")
            pdf(file = outputpath, width = 7,   height = 7)
            par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 
                                                                                                                
            color  = densCols(tmp_flu_corrected[tmp_gene], tmp_microevo_corrected_Kasasa[tmp_gene], colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
            plot(   tmp_flu_corrected[tmp_gene], tmp_microevo_corrected_Kasasa[tmp_gene], 
                    xlim    = c(-0.2, 0.8), 
                    ylim    = c(-0.2, 1.2), 
                    #xaxt   = "n",
                    #yaxt   = "n",
                    #log        = "xy",
                    xlab    = "Fluctuation (RM correction)",
                    ylab    = "Microevolutionary diversity (RM correction, inbred-Kasasa), Tech error (RM correction)",
                    pch     = 20,
                    cex     = 0.2,  #0.08,
                    col     = color, #rgb(0/255, 191/255, 255/255, alpha=0.1),  #"dodgerblue"
                    cex.lab = 1,
                    main    = paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
                    las     =1
                )
            
            mtext(paste(nrow(na.omit(cbind(tmp_flu_corrected[tmp_gene], tmp_microevo_corrected_Kasasa[tmp_gene]))), " genes", sep=""), side = 1, line = 6, adj = 1)

            cortest_list    <- cor.test(tmp_flu_corrected[tmp_gene], tmp_microevo_corrected_Kasasa[tmp_gene], method="spearman")
            mtext(paste("flu-evo rho = ", signif(cortest_list$estimate, 2), " (P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
            
            dev.off()
        }

         if(1){
            tmp_gene    <- intersect(names(tmp_flu_corrected), names(tmp_microevo_corrected_Oura))

            outputpath  <- paste(outputdir_tmp2, "RMCorrection_fluctuation-microevo(inbred-Oura)_", stage, ".pdf", sep="")
            pdf(file = outputpath, width = 7,   height = 7)
            par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 
                                                                                                                
            color  = densCols(tmp_flu_corrected[tmp_gene], tmp_microevo_corrected_Oura[tmp_gene], colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
            plot(   tmp_flu_corrected[tmp_gene], tmp_microevo_corrected_Oura[tmp_gene], 
                    xlim    = c(-0.2, 0.8), 
                    ylim    = c(-0.2, 1.2), 
                    #xaxt   = "n",
                    #yaxt   = "n",
                    #log        = "xy",
                    xlab    = "Fluctuation (RM correction)",
                    ylab    = "Microevolutionary diversity (RM correction, inbred-Oura), Tech error (RM correction)",
                    pch     = 20,
                    cex     = 0.2,  #0.08,
                    col     = color, #rgb(0/255, 191/255, 255/255, alpha=0.1),  #"dodgerblue"
                    cex.lab = 1,
                    main    = paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
                    las     =1
                )
            
            mtext(paste(nrow(na.omit(cbind(tmp_flu_corrected[tmp_gene], tmp_microevo_corrected_Oura[tmp_gene]))), " genes", sep=""), side = 1, line = 6, adj = 1)

            cortest_list    <- cor.test(tmp_flu_corrected[tmp_gene], tmp_microevo_corrected_Oura[tmp_gene], method="spearman")
            mtext(paste("flu-evo rho = ", signif(cortest_list$estimate, 2), " (P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
            
            dev.off()
        }

        if(1){
            outputpath  <- paste(outputdir_tmp2, "RMCorrection_explevel-microevo(inbred-Kasasa)_", stage, ".pdf", sep="")
            pdf(file = outputpath, width = 7,   height = 7)
            par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 

            color  = densCols(tmp_explevel[names(tmp_microevo_corrected_Kasasa)], tmp_microevo_corrected_Kasasa, colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
            plot(   tmp_explevel[names(tmp_microevo_corrected_Kasasa)], tmp_microevo_corrected_Kasasa, 
                    xlim    = c(0, 5), 
                    ylim    = c(-0.2, 1.2), 
                    #xaxt   = "n",
                    #yaxt   = "n",
                    #log        = "xy",
                    xlab    = "Expression level",
                    ylab    = "Microevolutionary diversity (RM correction, inbred-Kasasa)",
                    pch     = 20,
                    cex     = 0.1, #0.08,
                    col     = color, #"dodgerblue",
                    cex.lab = 1,
                    main    = paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
                    las     =1
                )
            
            mtext(paste(nrow(na.omit(cbind(tmp_explevel[names(tmp_microevo_corrected_Kasasa)], tmp_microevo_corrected_Kasasa))), " genes", sep=""), side = 1, line = 6, adj = 1)

            cortest_list    <- cor.test(tmp_explevel[names(tmp_microevo_corrected_Kasasa)], tmp_microevo_corrected_Kasasa, method="spearman")
            mtext(paste("microevo(RM correction)-exp level rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
            
            dev.off()
        }

        if(1){
            outputpath  <- paste(outputdir_tmp2, "RMCorrection_explevel-microevo(inbred-Oura)_", stage, ".pdf", sep="")
            pdf(file = outputpath, width = 7,   height = 7)
            par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 

            color  = densCols(tmp_explevel[names(tmp_microevo_corrected_Oura)], tmp_microevo_corrected_Oura, colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
            plot(   tmp_explevel[names(tmp_microevo_corrected_Oura)], tmp_microevo_corrected_Oura, 
                    xlim    = c(0, 5), 
                    ylim    = c(-0.2, 1.2), 
                    #xaxt   = "n",
                    #yaxt   = "n",
                    #log        = "xy",
                    xlab    = "Expression level",
                    ylab    = "Microevolutionary diversity (RM correction, inbred-Oura)",
                    pch     = 20,
                    cex     = 0.1, #0.08,
                    col     = color, #"dodgerblue",
                    cex.lab = 1,
                    main    = paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
                    las     =1
                )
            
            mtext(paste(nrow(na.omit(cbind(tmp_explevel[names(tmp_microevo_corrected_Oura)], tmp_microevo_corrected_Oura))), " genes", sep=""), side = 1, line = 6, adj = 1)

            cortest_list    <- cor.test(tmp_explevel[names(tmp_microevo_corrected_Oura)], tmp_microevo_corrected_Oura, method="spearman")
            mtext(paste("microevo(RM correction)-exp level rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
            
            dev.off()
        }

	}
}	



##########################################################################################################################
####　Distribution of developmental stability (no correction for expression level dependency) and gene expression levels			
##########################################################################################################################

outputdir_tmp2	<- paste(outputdir_tmp, "fluctuation-expressionlevel_distribution/", sep="")
dir.create(outputdir_tmp2)


# Developmental stability
if(1){
	fluctuation_early	<- na.omit(flu_nocorrection_early)
	fluctuation_mid		<- na.omit(flu_nocorrection_mid)
	fluctuation_late1	<- na.omit(flu_nocorrection_late1)
	fluctuation_late2	<- na.omit(flu_nocorrection_late2)


	## Plot only genes within the interquartile range
	get_data_inIQR	<- function(fluctuation_array){
		tmp	<- fluctuation_array[!(fluctuation_array %in% boxplot(fluctuation_array)$out)]
		return(tmp)
	}
	tmp_data    <- list(get_data_inIQR(fluctuation_early), get_data_inIQR(fluctuation_mid), get_data_inIQR(fluctuation_late1), get_data_inIQR(fluctuation_late2))


	## Obtain median value
	tmp_median		<- matrix(nrow=4,ncol=2)
	tmp_median[,1]	<- c(median(fluctuation_early), median(fluctuation_mid), median(fluctuation_late1), median(fluctuation_late2))
	tmp_median[,2]	<- c(1:4)
	
	output_path_plot			<- paste(outputdir_tmp2, "NORMcorrectedFluctuation_dotplot.pdf", sep = "")
	pdf(file = output_path_plot,	width = 9,		height = 9)
	par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(4, 1, 0), pty="s", yaxs="i")

	vioplot(	tmp_data,
				horizontal	= TRUE,
				xlim	= c(0.5,4.5),
				ylim	= c(min(unlist(tmp_data), na.rm=T),max(unlist(tmp_data), na.rm=T)),
				xlab	= "Fluctuation (NO RM correction)",
				ylab	= "", 
				names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
				rectCol	= NA,
				lineCol	= NA,
				border	= "grey35",
				colMed	= NA,
				main	= paste("After gene selection, Log, TPM=>",TPMCutOff, sep=""),
				las		= 1
				)

	par(new=T)

	plot(	tmp_median, 
			pch		= 16,
			col		= "gray60",
			xlim	= c(min(unlist(tmp_data), na.rm=T),max(unlist(tmp_data), na.rm=T)),
			ylim	= c(0.5,4.5),
			cex		= 2,
			xlab	= "",
			ylab	= "", 
			xaxt	= "n",
			yaxt	= "n"
		)

	axis(	side		= 4, 			
			at			= 1:4,  
    		labels		= c(paste("25%ile\t", signif(quantile(fluctuation_early,na.rm="TRUE")[2], 2), "\nmed\t", signif(quantile(fluctuation_early,na.rm="TRUE")[3], 2), "\n75%ile\t", signif(quantile(fluctuation_early,na.rm="TRUE")[4], 2), sep=""),
    						paste("25%ile\t", signif(quantile(fluctuation_mid,na.rm="TRUE")[2], 2), "\nmed\t", signif(quantile(fluctuation_mid,na.rm="TRUE")[3], 2), "\n75%ile\t", signif(quantile(fluctuation_mid,na.rm="TRUE")[4], 2), sep=""),
    						paste("25%ile\t", signif(quantile(fluctuation_late1,na.rm="TRUE")[2], 2), "\nmed\t", signif(quantile(fluctuation_late1,na.rm="TRUE")[3], 2), "\n75%ile\t", signif(quantile(fluctuation_late1,na.rm="TRUE")[4], 2), sep=""),
    						paste("25%ile\t", signif(quantile(fluctuation_late2,na.rm="TRUE")[2], 2), "\nmed\t", signif(quantile(fluctuation_late2,na.rm="TRUE")[3], 2), "\n75%ile\t", signif(quantile(fluctuation_late2,na.rm="TRUE")[4], 2), sep="")
    						),
    		tck			= 0,			
    		las			= 1,  			
   			mgp			= c(3.5,0.7,0), 	
   			cex.axis	= 0.8
		)

	# Kruskal-Walis test
	result <- kruskal.test(x=list(fluctuation_early,fluctuation_mid,fluctuation_late1,fluctuation_late2))
	mtext(paste("P = ", signif(result$p.value, 2), " (Kruskal-Walis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


	# Statistical test (Steel Dwass)
	dist				<- c(	fluctuation_early,fluctuation_mid,fluctuation_late1,fluctuation_late2)
	stage				<- c(	rep("Early",length(fluctuation_early)),
								rep("Mid",length(fluctuation_mid)),
								rep("Late1",length(fluctuation_late1)),
								rep("Late2",length(fluctuation_late2))
							)
	stage 				<- as.factor(stage)
	SD_pvalues			<- pSDCFlig(dist, stage, method = "Asymptotic")
	mtext(paste("Early-Mid(Steel-Dwass)", SD_pvalues$p.val[3],  sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
	mtext(paste("Early-Late1(Steel-Dwass)", SD_pvalues$p.val[1], sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
	mtext(paste("Early-Late2(Steel-Dwass)", SD_pvalues$p.val[2], sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
	mtext(paste("Mid-Late1(Steel-Dwass)", SD_pvalues$p.val[5], sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
	mtext(paste("Mid-Late2(Steel-Dwass)", SD_pvalues$p.val[6],  sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
	mtext(paste("Late1-Late2(Steel-Dwass)", SD_pvalues$p.val[4], sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)
						
	dev.off()

	rm(fluctuation_early)
	rm(fluctuation_mid)
	rm(fluctuation_late1)
	rm(fluctuation_late2)
	gc()
	gc()
}


# Expression levels
if(0){
	mean_early	<- apply(na.omit(TPM_early), 1, mean)
	mean_early	<- mean_early[names(flu_RMcorrection_early)]
	mean_mid	<- apply(na.omit(TPM_mid), 1, mean)
	mean_mid	<- mean_mid[names(flu_RMcorrection_mid)]
	mean_late1	<- apply(na.omit(TPM_late1), 1, mean)
	mean_late1	<- mean_late1[names(flu_RMcorrection_late1)]
	mean_late2	<- apply(na.omit(TPM_late2), 1, mean)
	mean_late2	<- mean_late2[names(flu_RMcorrection_late2)]

	
	output_path_plot			<- paste(outputdir_tmp2, "expressionlevel_dotplot.pdf", sep = "")
	pdf(file = output_path_plot,	width = 9,		height = 9)
	par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(4, 1, 0), pty="s")

	beeswarm(	list(mean_early,mean_mid,mean_late1,mean_late2),
				col 	= c(bg=rgb(0,154/255,205/255, alpha=0.9),
 							bg=rgb(151/255,255/255,255/255, alpha=0.9),
 							bg=rgb(144/255,238/255,144/255, alpha=0.9), 
							bg=rgb(255/255,193/255,37/255,  alpha=0.9)
							),
				log		= "TRUE", 
				xlim	= c(0.0001,5),
				spacing = 0.6,
				xaxt	= "n",
				yaxt	= "n",
				vertical = F,
				pch		=16, 
				cex		= 0.2
			)

	par(new=T)

	boxplot(	mean_early,mean_mid,mean_late1,mean_late2,
				border	= "grey28", #"grey35",
				names	= c("Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"	),
				horizontal	= T,
				boxwex	= 0.5,
				col		= NA,
				log		= "x",
				ylim	= c(0.001,5),
				xlab	= "Expression level",
				cex.lab	= 1,
				cex.axis	= 1,
				cex 	= 0.1,
				pch 	= 16,
				las		= 1,
				main	= paste("After gene selection, Log, meanTPM=>",TPMCutOff, sep="")
			)


	axis(	side		= 4, 	
			at			= 1:4,  
    		labels		= c(paste("25%ile\t", signif(quantile(mean_early,na.rm="TRUE")[2], 2), "\nmed\t", signif(quantile(mean_early,na.rm="TRUE")[3], 2), "\n75%ile\t", signif(quantile(mean_early,na.rm="TRUE")[4], 2), sep=""),
    						paste("25%ile\t", signif(quantile(mean_mid,na.rm="TRUE")[2], 2), "\nmed\t", signif(quantile(mean_mid,na.rm="TRUE")[3], 2), "\n75%ile\t", signif(quantile(mean_mid,na.rm="TRUE")[4], 2), sep=""),
    						paste("25%ile\t", signif(quantile(mean_late1,na.rm="TRUE")[2], 2), "\nmed\t", signif(quantile(mean_late1,na.rm="TRUE")[3], 2), "\n75%ile\t", signif(quantile(mean_late1,na.rm="TRUE")[4], 2), sep=""),
    						paste("25%ile\t", signif(quantile(mean_late2,na.rm="TRUE")[2], 2), "\nmed\t", signif(quantile(mean_late2,na.rm="TRUE")[3], 2), "\n75%ile\t", signif(quantile(mean_late2,na.rm="TRUE")[4], 2), sep="")
    						),
    		tck			= 0,			
    		las			= 1,  		
   			mgp			= c(3.5,0.7,0), 	
   			cex.axis	= 0.8
		)



	# Kruskal-Walis test
	result <- kruskal.test(x=list(mean_early,mean_mid,mean_late1,mean_late2))
	mtext(paste("P = ", signif(result$p.value, 2), " (Kruskal-Walis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)
	

	# Statistical test (Steel Dwass)
	dist				<- c(	mean_early,mean_mid,mean_late1,mean_late2)
	stage				<- c(	rep("Early",length(mean_early)),
									rep("Mid",length(mean_mid)),
									rep("Late1",length(mean_late1)),
									rep("Late2",length(mean_late2))
							)
	stage 				<- as.factor(stage)
	SD_pvalues			<- pSDCFlig(dist, stage, method = "Asymptotic")
	mtext(paste("Early-Mid(Steel-Dwass)", SD_pvalues$p.val[3],  sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
	mtext(paste("Early-Late1(Steel-Dwass)", SD_pvalues$p.val[1], sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
	mtext(paste("Early-Late2(Steel-Dwass)", SD_pvalues$p.val[2], sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
	mtext(paste("Mid-Late1(Steel-Dwass)", SD_pvalues$p.val[5], sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
	mtext(paste("Mid-Late2(Steel-Dwass)", SD_pvalues$p.val[6],  sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
	mtext(paste("Late1-Late2(Steel-Dwass)", SD_pvalues$p.val[4], sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)
						
	dev.off()


	rm(mean_early)
	rm(mean_mid)
	rm(mean_late1)
	rm(mean_late2)
	gc()
	gc()
}




##############################################################################################################
#### Pleiotropy analysis
##############################################################################################################

outputdir_tmp2	<- paste(outputdir_tmp, "Pleiotropy/", sep="")
dir.create(outputdir_tmp2)

# Calculate the number of stages and tissues that express the gene for each gene.
if(1){
	## Tissue
	if(1){
		### Read gene expression data by tissue
		#num_tissue_EachGene_expressing		<- read.table("~/Documents/work/Research/Fluctuation/Data/ExpressionTables/num_tissue_EachGene_expressing.txt", header=TRUE, sep="\t", row.names=1)
		TissueTPM_tmp 			<-  read.table("/Volumes/HDD_BackUp/Fluctuation/Data/Medaka_AdultTissue_TPM.txt", header=TRUE, sep="\t", row.names=1)

		### Calculate average expression among 4 replicates
		TPM_tissue_mean				<- matrix(nrow = nrow(TissueTPM_tmp), ncol = 25)
		colnames(TPM_tissue_mean)	<- rep(NA,25)
		rownames(TPM_tissue_mean)	<- rownames(TissueTPM_tmp)
		for (i in 1:25){ 
			TPM_tissue_mean[,i]	<- apply(TissueTPM_tmp[,(4*i-3):(4*i)], 1, mean)
			colnames(TPM_tissue_mean)[i]	<- strsplit(colnames(TissueTPM_tmp)[(4*i-3)],"_")[[1]][1]
		}
		TPM_tissue_mean			<- log(TPM_tissue_mean+1, base=10)


		### Converted to expressed or not expressed
		cutoff_function	<- function(x){
			if (x >= TPMCutOff_TissueStage){
				if(x > 0){
					return (1)
				}else{
					return (0)
				}
				
			}else{
				return (0)
			}
		}
		TPM_tissue_mean	<- apply(TPM_tissue_mean, c(1,2), cutoff_function)

	
		### Calculate the number of tissues expressed
		num_tissue_EachGene_expressing <- data.frame(apply(TPM_tissue_mean, 1, function(x) sum(x, na.rm=TRUE)))

		rm(TissueTPM_tmp)
		rm(TPM_tissue_mean)
		gc()
		gc()
	}

	## Stage
	if(1){
		### Read time-course gene expression data
		num_stage_EachGene_expressing_tmp	<- read.table("~/Documents/work/Research/Fluctuation/Data/ExpressionTables/temporal_pleiotropy/Medaka_TPM_TimeCourse.txt", header=TRUE, sep="\t", row.names=1)

		TPM_timecourse_mean					<- matrix(nrow = nrow(num_stage_EachGene_expressing_tmp), ncol = 16)
		for (i in 1:16){ 
			TPM_timecourse_mean[,i]			<- apply(num_stage_EachGene_expressing_tmp[,(3*i-2):(3*i)], 1, mean)
		}
		TPM_timecourse_mean					<- log(TPM_timecourse_mean+1, base=10)

		convert_1or0	<- function(x){
			if(x >= TPMCutOff_TissueStage){
				if(x > 0){
					return (1)
				}else{
					return (0)
				}
			}else{
				return (0)
			}
		}

		TPM_1or0_temporal							<- apply(TPM_timecourse_mean, c(1,2), convert_1or0)
		num_stages_EachGene_expressing				<- apply(TPM_1or0_temporal, 1, sum)
		num_stages_EachGene_expressing 				<- data.frame(num_stages_EachGene_expressing)
		rownames(num_stages_EachGene_expressing)	<- rownames(num_stage_EachGene_expressing_tmp)

		#outputpath	<- paste(outputdir_tmp2, "num_stages_EachGene_expressing.txt", sep="")
		#write.table(num_stages_EachGene_expressing, file = outputpath, quote = FALSE, sep = "\t", col.names = FALSE, row.names = TRUE )

		rm(num_stage_EachGene_expressing_tmp)
		gc()
		gc()
	}

}

## Visualization
if(1){
	stages	<- c("early", "mid", "late1", "late2")
	for(i in 1:4){
		stage		<- stages[i]
		tmp_inbred	<- na.omit(get(paste("TPM_", stage, sep="")))

		tmp_flu_corrected						<- get(paste("flu_RMcorrection_", stage, sep=""))

		tmp_microevo_corrected_new				<- get(paste("microevo_RMcorrection_new_", stage, sep=""))

		tmp_genes								<- intersect(rownames(tmp_inbred), names(tmp_flu_corrected))
		tmp_explevel							<- apply(tmp_inbred[tmp_genes,], 1, mean)	


		# Function to get genes within the interquartile range
		get_data_inIQR	<- function(k){
			if(length(tmp_box_table[tmp_box_table[,2] == k,1]) == 0){
				return (NA)
			}else{
				return(tmp_box_table[tmp_box_table[,2] == k,1][!(tmp_box_table[tmp_box_table[,2] == k,1] %in% boxplot(tmp_box_table[tmp_box_table[,2] == k,1])$out)])
			}
		}
		
		get_wex			<- function(k){
			if(length(tmp_box_table[tmp_box_table[,2] == k,1]) == 0){
				return ()
			}else{
				return(length(tmp_box_table[tmp_box_table[,2] == k,1][!(tmp_box_table[tmp_box_table[,2] == k,1] %in% boxplot(tmp_box_table[tmp_box_table[,2] == k,1])$out)]))
			}
		}


		## Plot (developmental stability vs spatial pleiotropy)
		if(1){
			tmp_box_table		<- cbind(tmp_flu_corrected, num_tissue_EachGene_expressing[names(tmp_flu_corrected),1])


            tmp_data    <- list(get_data_inIQR(1), get_data_inIQR(2), get_data_inIQR(3), get_data_inIQR(4), get_data_inIQR(5),
								get_data_inIQR(6), get_data_inIQR(7), get_data_inIQR(8), get_data_inIQR(9), get_data_inIQR(10), 
								get_data_inIQR(11), get_data_inIQR(12), get_data_inIQR(13), get_data_inIQR(14), get_data_inIQR(15),
								get_data_inIQR(16), get_data_inIQR(17), get_data_inIQR(18), get_data_inIQR(19), get_data_inIQR(20),
								get_data_inIQR(21), get_data_inIQR(22), get_data_inIQR(23),  get_data_inIQR(24), get_data_inIQR(25)
                            )

			tmp_wex		<- c(	get_wex(1), get_wex(2), get_wex(3), get_wex(4), get_wex(5),
								get_wex(6), get_wex(7), get_wex(8), get_wex(9), get_wex(10), 
								get_wex(11), get_wex(12), get_wex(13), get_wex(14), get_wex(15),
								get_wex(16), get_wex(17), get_wex(18), get_wex(19), get_wex(20),
								get_wex(21), get_wex(22), get_wex(23), get_wex(24), get_wex(25)
                            )
			tmp_wex		<- (tmp_wex/sum(tmp_wex))*2

			tmp_median		<- matrix(nrow=25,ncol=2)
			tmp_median[,1]	<- c(1:25)
			tmp_median[,2]	<- c(	median(tmp_box_table[tmp_box_table[,2] == 1,1]),
									median(tmp_box_table[tmp_box_table[,2] == 2,1]),
									median(tmp_box_table[tmp_box_table[,2] == 3,1]),
									median(tmp_box_table[tmp_box_table[,2] == 4,1]),
									median(tmp_box_table[tmp_box_table[,2] == 5,1]),
									median(tmp_box_table[tmp_box_table[,2] == 6,1]),
									median(tmp_box_table[tmp_box_table[,2] == 7,1]),
									median(tmp_box_table[tmp_box_table[,2] == 8,1]),
									median(tmp_box_table[tmp_box_table[,2] == 9,1]),
									median(tmp_box_table[tmp_box_table[,2] == 10,1]),
									median(tmp_box_table[tmp_box_table[,2] == 11,1]),
									median(tmp_box_table[tmp_box_table[,2] == 12,1]),
									median(tmp_box_table[tmp_box_table[,2] == 13,1]),
									median(tmp_box_table[tmp_box_table[,2] == 14,1]),
									median(tmp_box_table[tmp_box_table[,2] == 15,1]),
									median(tmp_box_table[tmp_box_table[,2] == 16,1]),
									median(tmp_box_table[tmp_box_table[,2] == 17,1]),
									median(tmp_box_table[tmp_box_table[,2] == 18,1]),
									median(tmp_box_table[tmp_box_table[,2] == 19,1]),
									median(tmp_box_table[tmp_box_table[,2] == 20,1]),
									median(tmp_box_table[tmp_box_table[,2] == 21,1]),
									median(tmp_box_table[tmp_box_table[,2] == 22,1]),
									median(tmp_box_table[tmp_box_table[,2] == 23,1]),
									median(tmp_box_table[tmp_box_table[,2] == 24,1]),
									median(tmp_box_table[tmp_box_table[,2] == 25,1])
								)

			plot_position	<- c(1:25)
			plot_position	<- plot_position[unlist(lapply(tmp_data, function(x) !(is.na(sum(x, na.rm=F)))))]


			outputpath	<- paste(outputdir_tmp2, "RMCorrection_fluctuation-tissue_", stage, ".pdf", sep="")
			pdf(file = outputpath, width = 7,	height = 7)
			par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", yaxs="i") # 

			tmp_data	<- lapply(tmp_data, na.omit)
			vioplot(	tmp_data,			
						main	= paste("After gene selection, Log, TPM=>",TPMCutOff, ", (outlitter = F), ", stage, sep=""),
						#areaEqual = T, 
						xlim	= c(0.5,25.5),
						ylim	= c(min(unlist(tmp_data), na.rm=T),max(unlist(tmp_data), na.rm=T)),
						at		= plot_position,
						xaxt	= "n",
						xlab	= "Number of tissues",
						ylab	= "Fluctuation (RM correction)", 
						col		= "turquoise2", 
						rectCol	= NA,
						lineCol	= NA,
						border	= "turquoise4",
						colMed	= NA,
						wex		= tmp_wex,
						las		= 1
					)

			axis(1, at=c(1:25), labels=c(1,rep(NA,3), 5, rep(NA,4), 10, rep(NA,4), 15, rep(NA,4), 20, rep(NA,4), 25))

			par(new=T)

			plot(	tmp_median, 
					pch		= 16,
					col		= "turquoise4",
					xlim	= c(0.5,25.5),
					ylim	= c(min(unlist(tmp_data), na.rm=T),max(unlist(tmp_data), na.rm=T)),
					xlab	= "",
					ylab	= "", 
					xaxt	= "n",
					yaxt	= "n"
					)

				cortest_list	<- cor.test(tmp_box_table[,1], tmp_box_table[,2], method = "spearman")
				mtext(paste("fluctuation(RM correction)-tissue rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 6, adj = 1)

				# Statistical test Kruskal Wallis)
				x=list(	if(length(tmp_box_table[tmp_box_table[,2] == 1,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 1,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 2,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 2,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 3,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 3,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 4,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 4,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 5,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 5,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 6,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 6,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 7,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 7,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 8,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 8,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 9,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 9,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 10,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 10,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 11,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 11,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 12,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 12,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 13,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 13,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 14,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 14,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 15,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 15,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 16,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 16,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 17,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 17,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 18,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 18,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 19,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 19,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 20,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 20,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 21,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 21,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 22,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 22,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 23,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 23,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 24,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 24,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 25,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 25,1]}
				)
				KW_result <- kruskal.test(x[!is.na(x)])
				mtext(paste("p-values (Kruskal-Wallis)    ", KW_result$p.value, sep=""), side = 1, line = 7, adj = 1)
			
				dev.off()
		}

		## Plot (developmental stability vs temporal pleiotropy)
		if(1){
			tmp_box_table		<- merge(tmp_flu_corrected, num_stages_EachGene_expressing,by="row.names")
			tmp_box_table		<- tmp_box_table[,c(2,3)]

            tmp_data    <- list(get_data_inIQR(1), get_data_inIQR(2), get_data_inIQR(3), get_data_inIQR(4),
								get_data_inIQR(5), get_data_inIQR(6), get_data_inIQR(7), get_data_inIQR(8),
								get_data_inIQR(9), get_data_inIQR(10), get_data_inIQR(11), get_data_inIQR(12), 
								get_data_inIQR(13), get_data_inIQR(14), get_data_inIQR(15),	get_data_inIQR(16)
                            )

			tmp_wex		<- c(	get_wex(1), get_wex(2), get_wex(3), get_wex(4),
								get_wex(5), get_wex(6), get_wex(7), get_wex(8),
								get_wex(9), get_wex(10), get_wex(11), get_wex(12), 
								get_wex(13), get_wex(14), get_wex(15),	get_wex(16)
                            )

            tmp_wex     <- (tmp_wex/sum(tmp_wex))*2

            tmp_median      <- matrix(nrow=16,ncol=2)
            tmp_median[,1]  <- c(1:16)
            tmp_median[,2]  <- c(   median(tmp_box_table[tmp_box_table[,2] == 1,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 2,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 3,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 4,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 5,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 6,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 7,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 8,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 9,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 10,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 11,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 12,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 13,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 14,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 15,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 16,1])
                                )

			plot_position	<- c(1:16)
			plot_position	<- plot_position[unlist(lapply(tmp_data, function(x) !(is.na(sum(x, na.rm=F)))))]


			outputpath	<- paste(outputdir_tmp2, "RMCorrection_fluctuation-stage_", stage, ".pdf", sep="")
			pdf(file = outputpath, width = 7,	height = 7)
			par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", yaxs="i") 

			tmp_data	<- lapply(tmp_data, na.omit)
            vioplot(    tmp_data,           
                        main    = paste("After gene selection, Log, TPM=>",TPMCutOff, ", (outlitter = F), ", stage, sep=""),
                        #areaEqual = T, 
                        xlim    = c(0.5,16.5),
                        ylim    = c(min(unlist(tmp_data), na.rm=T),max(unlist(tmp_data), na.rm=T)),
						at		= plot_position,
						xaxt	= "n",
                        xlab    = "Number of stages",
                        ylab    = "Fluctuation (RM correction)", 
                        col     = "indianred2", 
                        rectCol = NA,
                        lineCol = NA,
                        border  = "indianred4",
                        colMed  = NA,
                        wex     = tmp_wex,
                        las     = 1
                    )

			axis(1, at=c(1:16), labels=c(1,rep(NA,2), 4, rep(NA,3), 8, rep(NA,3), 12, rep(NA,3), 16))

            par(new=T)

            plot(   tmp_median, 
                    pch     = 16,
                    col     = "indianred4",
                    xlim    = c(0.5,16.5),
                    ylim    = c(min(unlist(tmp_data), na.rm=T),max(unlist(tmp_data), na.rm=T)),
                    xlab    = "",
                    ylab    = "", 
                    xaxt    = "n",
                    yaxt    = "n"
                    )

			cortest_list	<- cor.test(tmp_box_table[,1], tmp_box_table[,2], method = "spearman")
			mtext(paste("fluctuation-stage rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 6, adj = 1)

			# Statistical test Kruskal Wallis)
			x=list(	if(length(tmp_box_table[tmp_box_table[,2] == 1,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 1,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 2,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 2,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 3,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 3,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 4,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 4,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 5,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 5,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 6,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 6,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 7,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 7,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 8,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 8,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 9,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 9,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 10,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 10,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 11,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 11,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 12,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 12,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 13,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 13,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 14,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 14,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 15,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 15,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 16,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 16,1]}
				)
			KW_result <- kruskal.test(x[!is.na(x)])
			mtext(paste("p-values (Kruskal-Wallis)    ", KW_result$p.value, sep=""), side = 1, line = 7, adj = 1)
			
			dev.off()
		}



		## Plot (intraspecies diversity vs spatial pleiotropy)
		if(1){
			tmp_box_table		<- cbind(tmp_microevo_corrected_new, num_tissue_EachGene_expressing[names(tmp_microevo_corrected_new),1])

            tmp_data    <- list(get_data_inIQR(1), get_data_inIQR(2), get_data_inIQR(3), get_data_inIQR(4), get_data_inIQR(5),
								get_data_inIQR(6), get_data_inIQR(7), get_data_inIQR(8), get_data_inIQR(9), get_data_inIQR(10), 
								get_data_inIQR(11), get_data_inIQR(12), get_data_inIQR(13), get_data_inIQR(14), get_data_inIQR(15),
								get_data_inIQR(16), get_data_inIQR(17), get_data_inIQR(18), get_data_inIQR(19), get_data_inIQR(20),
								get_data_inIQR(21), get_data_inIQR(22), get_data_inIQR(23),  get_data_inIQR(24), get_data_inIQR(25)
                            )

			tmp_wex		<- c(	get_wex(1), get_wex(2), get_wex(3), get_wex(4), get_wex(5),
								get_wex(6), get_wex(7), get_wex(8), get_wex(9), get_wex(10), 
								get_wex(11), get_wex(12), get_wex(13), get_wex(14), get_wex(15),
								get_wex(16), get_wex(17), get_wex(18), get_wex(19), get_wex(20),
								get_wex(21), get_wex(22), get_wex(23), get_wex(24), get_wex(25)
							)	
            tmp_wex     <- (tmp_wex/sum(tmp_wex))*2

            tmp_median      <- matrix(nrow=25,ncol=2)
            tmp_median[,1]  <- c(1:25)
            tmp_median[,2]  <- c(   median(tmp_box_table[tmp_box_table[,2] == 1,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 2,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 3,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 4,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 5,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 6,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 7,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 8,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 9,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 10,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 11,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 12,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 13,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 14,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 15,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 16,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 17,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 18,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 19,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 20,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 21,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 22,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 23,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 24,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 25,1])
                                )

			plot_position	<- c(1:25)
			plot_position	<- plot_position[unlist(lapply(tmp_data, function(x) !(is.na(sum(x, na.rm=F)))))]


			outputpath	<- paste(outputdir_tmp2, "RMCorrection_microevo(new)-tissue_", stage, ".pdf", sep="")
			pdf(file = outputpath, width = 7,	height = 7)
			par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", yaxs="i") # 

			tmp_data	<- lapply(tmp_data, na.omit)
            vioplot(    tmp_data,           
                        main    = paste("After gene selection, Log, TPM=>",TPMCutOff, ", (outlitter = F), ", stage, sep=""),
                        #areaEqual = T, 
                        xlim    = c(0.5,25.5),
                        ylim    = c(min(unlist(tmp_data), na.rm=T),max(unlist(tmp_data), na.rm=T)),
						at		= plot_position,
						xaxt	= "n",
                        xlab    = "Number of tissues",
                        ylab    = "Microevolutionary diversity (RM correction, new)", 
                        col     = "turquoise2", 
                        rectCol = NA,
                        lineCol = NA,
                        border  = "turquoise4",
                        colMed  = NA,
                        wex     = tmp_wex,
                        las     = 1
                    )

			axis(1, at=c(1:25), labels=c(1,rep(NA,3), 5, rep(NA,4), 10, rep(NA,4), 15, rep(NA,4), 20, rep(NA,4), 25))

            par(new=T)

            plot(   tmp_median, 
                    pch     = 16,
                    col     = "turquoise4",
                    xlim    = c(0.5,25.5),
                    ylim    = c(min(unlist(tmp_data), na.rm=T),max(unlist(tmp_data), na.rm=T)),
                    xlab    = "",
                    ylab    = "", 
                    xaxt    = "n",
                    yaxt    = "n"
                    )

			cortest_list	<- cor.test(tmp_box_table[,1], tmp_box_table[,2], method = "spearman")
			mtext(paste("microevo(new)-tissue rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 6, adj = 1)

			# Statistical test Kruskal Wallis)
			x=list(	if(length(tmp_box_table[tmp_box_table[,2] == 1,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 1,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 2,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 2,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 3,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 3,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 4,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 4,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 5,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 5,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 6,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 6,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 7,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 7,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 8,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 8,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 9,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 9,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 10,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 10,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 11,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 11,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 12,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 12,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 13,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 13,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 14,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 14,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 15,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 15,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 16,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 16,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 17,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 17,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 18,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 18,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 19,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 19,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 20,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 20,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 21,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 21,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 22,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 22,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 23,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 23,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 24,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 24,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 25,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 25,1]}
				)
			KW_result <- kruskal.test(x[!is.na(x)])
			mtext(paste("p-values (Kruskal-Wallis)    ", KW_result$p.value, sep=""), side = 1, line = 7, adj = 1)
			
			dev.off()
		}



		## Plot (intraspecies diversity vs temporal pleiotropy)
		if(1){
			tmp_box_table		<- merge(tmp_microevo_corrected_new, num_stages_EachGene_expressing, by="row.names")
			tmp_box_table		<- tmp_box_table[,c(2,3)]

            tmp_data    <- list(get_data_inIQR(1), get_data_inIQR(2), get_data_inIQR(3), get_data_inIQR(4),
								get_data_inIQR(5), get_data_inIQR(6), get_data_inIQR(7), get_data_inIQR(8),
								get_data_inIQR(9), get_data_inIQR(10), get_data_inIQR(11), get_data_inIQR(12), 
								get_data_inIQR(13), get_data_inIQR(14), get_data_inIQR(15),	get_data_inIQR(16)
                            )

			tmp_wex		<- c(	get_wex(1), get_wex(2), get_wex(3), get_wex(4),
								get_wex(5), get_wex(6), get_wex(7), get_wex(8),
								get_wex(9), get_wex(10), get_wex(11), get_wex(12), 
								get_wex(13), get_wex(14), get_wex(15),	get_wex(16)
                            )

            tmp_wex     <- (tmp_wex/sum(tmp_wex))*2

            tmp_median      <- matrix(nrow=16,ncol=2)
            tmp_median[,1]  <- c(1:16)
            tmp_median[,2]  <- c(   median(tmp_box_table[tmp_box_table[,2] == 1,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 2,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 3,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 4,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 5,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 6,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 7,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 8,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 9,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 10,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 11,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 12,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 13,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 14,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 15,1]),
                                    median(tmp_box_table[tmp_box_table[,2] == 16,1])
                                )

			plot_position	<- c(1:16)
			plot_position	<- plot_position[unlist(lapply(tmp_data, function(x) !(is.na(sum(x, na.rm=F)))))]


			outputpath	<- paste(outputdir_tmp2, "RMCorrection_microevo(new)-stage_", stage, ".pdf", sep="")
			pdf(file = outputpath, width = 7,	height = 7)
			par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", yaxs="i") # 

			tmp_data	<- lapply(tmp_data, na.omit)
            vioplot(    tmp_data,           
                        main    = paste("After gene selection, Log, TPM=>",TPMCutOff, ", (outlitter = F), ", stage, sep=""),
                        #areaEqual = T, 
                        xlim    = c(0.5,16.5),
                        ylim    = c(min(unlist(tmp_data), na.rm=T),max(unlist(tmp_data), na.rm=T)),
						at		= plot_position,
						xaxt	= "n",
                        xlab    = "Number of stages",
                        ylab    = "Microevolutionary diversity (RM correction, new)", 
                        col     = "indianred2", 
                        rectCol = NA,
                        lineCol = NA,
                        border  = "indianred4",
                        colMed  = NA,
                        wex     = tmp_wex,
                        las     = 1
                    )

			axis(1, at=c(1:16), labels=c(1,rep(NA,2), 4, rep(NA,3), 8, rep(NA,3), 12, rep(NA,3), 16))

            par(new=T)

            plot(   tmp_median, 
                    pch     = 16,
                    col     = "indianred4",
                    xlim    = c(0.5,16.5),
                    ylim    = c(min(unlist(tmp_data)),max(unlist(tmp_data))),
                    xlab    = "",
                    ylab    = "", 
                    xaxt    = "n",
                    yaxt    = "n"
                    )

			cortest_list	<- cor.test(tmp_box_table[,1], tmp_box_table[,2], method = "spearman")
			mtext(paste("microevo(new)-stage rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 6, adj = 1)

			# Statistical test Kruskal Wallis)
			x=list(	if(length(tmp_box_table[tmp_box_table[,2] == 1,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 1,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 2,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 2,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 3,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 3,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 4,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 4,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 5,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 5,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 6,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 6,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 7,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 7,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 8,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 8,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 9,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 9,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 10,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 10,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 11,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 11,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 12,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 12,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 13,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 13,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 14,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 14,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 15,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 15,1]},
					if(length(tmp_box_table[tmp_box_table[,2] == 16,1]) == 0){NA}else{tmp_box_table[tmp_box_table[,2] == 16,1]}
				)
			KW_result <- kruskal.test(x[!is.na(x)])
			mtext(paste("p-values (Kruskal-Wallis)    ", KW_result$p.value, sep=""), side = 1, line = 7, adj = 1)
			
			dev.off()
		}


		tmp				<- cbind(	tmp_flu_corrected, 
									tmp_microevo_corrected_new[names(tmp_flu_corrected)], 
									num_tissue_EachGene_expressing[names(tmp_flu_corrected),1], 
									tmp_explevel[names(tmp_flu_corrected)]
									)
		colnames(tmp)	<- c("fluctuation", "microevo_diversity_new", "spatial_pleiotropy", "exp_level")
		cor_tmp			<- pcor(na.omit(tmp), method="spearman")
		cor_tmp			<- rbind(cor_tmp$estimate, c("P-value","","",""), cor_tmp$p.value)
		outputpath		<- paste(outputdir_tmp2, "all(new)_tissue_partialcor_", stage, ".txt", sep="")
		write.table(cor_tmp, file = outputpath, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE )


		tmp_gene		<- intersect(rownames(num_stages_EachGene_expressing), names(tmp_flu_corrected))
		tmp				<- cbind(	tmp_flu_corrected[tmp_gene], 
									tmp_microevo_corrected_new[tmp_gene],
									num_stages_EachGene_expressing[tmp_gene,],
									tmp_explevel[tmp_gene]
									)
		colnames(tmp)	<- c("fluctuation", "microevo_diversity_new",  "temporal_pleiotropy", "exp_level")
		cor_tmp			<- pcor(na.omit(tmp), method="spearman")
		cor_tmp			<- rbind(cor_tmp$estimate, c("P-value","","",""), cor_tmp$p.value)
		outputpath		<- paste(outputdir_tmp2, "all(new)_stage_partialcor_", stage, ".txt", sep="")
		write.table(cor_tmp, file = outputpath, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE )
	}
}




##############################################################################################################
#### Developmental stability vs interspecies diversity
##############################################################################################################

outputdir_tmp2	<- paste(outputdir_tmp, "Macroevo/", sep="")
dir.create(outputdir_tmp2)


if(1){
	tmp_inbred				<- na.omit(TPM_mid)

	tmp_flu_corrected						<- flu_RMcorrection_mid

	tmp_microevo_corrected_conventional		<- microevo_RMcorrection_conventional_mid
	tmp_microevo_corrected_new				<- microevo_RMcorrection_new_mid
}


# vs Zebrafish
if(1){
	## Calculate intraspecies diversity and correct for expression level dependency
	if(1){
		RBBH_table_zebra					<- read.table("~/Documents/work/Research/Fluctuation/Data/ExpressionTables/RBBH/RBBH_gene_Oryzias_latipes_and_Danio_rerio.txt", header=FALSE, sep="\t")
		TPM_zebra							<- read.table("/Volumes/HDD_BackUp/Fluctuation/Data/Dr_tpm.coding.txt", header=TRUE, sep="\t", row.names = 1)
		TPM_zebra							<- log(TPM_zebra+1, base=10)

		## Obtain reciprocal blast best hit table
		tmp_ortholog						<- data.frame(RBBH_table_zebra[,1])
		rownames(tmp_ortholog)				<- tmp_ortholog[,1]
		medaka_ortholog_TPM					<- merge(tmp_inbred, tmp_ortholog, by.x = "row.names", by.y = "row.names")
		rownames(medaka_ortholog_TPM)		<- medaka_ortholog_TPM[,1]
		medaka_ortholog_TPM					<- medaka_ortholog_TPM[,2:(ncol(medaka_ortholog_TPM)-1)]	

		TPM_zebra2							<- TPM_zebra[,grep("prim5", colnames(TPM_zebra))]
		tmp_ortholog						<- data.frame(RBBH_table_zebra)
		rownames(tmp_ortholog)				<- tmp_ortholog[,2]
		zebra_ortholog_TPM					<- merge(TPM_zebra2, tmp_ortholog, by.x = "row.names", by.y = "row.names")
		rownames(zebra_ortholog_TPM)		<- zebra_ortholog_TPM[,5]
		zebra_ortholog_TPM					<- zebra_ortholog_TPM[,2:4]	

		combined_RBBH						<- merge(medaka_ortholog_TPM, zebra_ortholog_TPM, by.x = "row.names", by.y = "row.names")
		rownames(combined_RBBH)				<- combined_RBBH[,1]
		combined_RBBH						<- combined_RBBH[,-1]
		combined_RBBH						<- na.omit(combined_RBBH)


		### Calculate interspecies diversity 
		GXP_change_nocorrection			<- apply(combined_RBBH, 1, get_GXPchange_log)
		GXP_change_nocorrection			<- data.frame(GXP_change_nocorrection)
		tmp_gene						<- intersect(names(tmp_flu_corrected), rownames(GXP_change_nocorrection))

		tmp_explevel					<- apply(tmp_inbred[tmp_gene,], 1, mean)
		tmp_explevel					<- tmp_explevel[tmp_gene]
		tmp_explevel					<- tmp_explevel[order(tmp_explevel, decreasing=TRUE)]
		GXP_change_nocorrection			<- GXP_change_nocorrection[names(tmp_explevel),]
		names(GXP_change_nocorrection)	<- names(tmp_explevel)

		RM_GXP_change			<- numeric(length(GXP_change_nocorrection))
		names(RM_GXP_change)	<- names(GXP_change_nocorrection)
		for(k in 1:length(GXP_change_nocorrection)){
	
			if(k == 1 || k == length(GXP_change_nocorrection)){
				RM_GXP_change[k]	<- GXP_change_nocorrection[k]
			}else{
			
				if(k <= 250){
					RM_GXP_change[k]	<- median(GXP_change_nocorrection[1:(2*k-1)], na.rm=TRUE)
				}else{
					if(k <= (length(GXP_change_nocorrection) - 250)){
						RM_GXP_change[k]	<- median(GXP_change_nocorrection[(k-250):(k+250)], na.rm=TRUE)
					}else{
						RM_GXP_change[k]	<- median(GXP_change_nocorrection[(2*k-length(GXP_change_nocorrection)):length(GXP_change_nocorrection)], na.rm=TRUE)
					}						
				}					
			}	
		}

		GXP_change_corrected	<- GXP_change_nocorrection - RM_GXP_change
		GXP_change_corrected	<- data.frame(GXP_change_corrected)
	}


	## Plot (developmental stability vs intraspecies diversity)
	if(1){
		merged				<- merge(tmp_flu_corrected, GXP_change_corrected, by.x = "row.names", by.y = "row.names") 
		rownames(merged)	<- merged[,1]
		merged				<- merged[,-1]

		outputpath	<- paste(outputdir_tmp2, "RMCorrection_fluctuation-macroevo_zebrafish.pdf", sep="")
		pdf(file = outputpath,	width = 7,		height = 7)
		par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i")
		#par(xpd = TRUE, mar = c(8, 6, 2, 2), mgp = c(5, 1, 0), pty="s")

		color  = densCols(merged[,1], merged[,2], colramp = colorRampPalette(c("#01579B","#0277BD","#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
		plot(	merged[,1], merged[,2],
				xlab		="Fluctuation (RM correction)", #Fluctuation
				ylab		= "Interspecies diversity (RM correction)", #Macroevolutionary change
				#log			= "xy",
				#xaxt	= "n",
				#yaxt	= "n",
				xlim		= c(-0.075, 0.3), 
				ylim		= c(-0.75, 2), 
				pch			= 16,
				col			= color, #"dodgerblue",
				cex			= 0.3,
				cex.lab		= 1,
				main		= paste("After gene selection, Log, TPM=>",TPMCutOff, ", zebrafish", sep=""),
				las			= 1
			)
		
		mtext(paste(nrow(na.omit(merged)), " genes", sep=""), side = 1, line = 6, adj = 1)
		cortest_list	<- cor.test(merged[,1], merged[,2], method = "spearman")
		mtext(paste("fluctuation-macroevo(RM correction for both) rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
	
		dev.off()																	
	}


	rm(RBBH_table_zebra)
	rm(TPM_zebra)
	rm(TPM_zebra2)
	rm(zebra_ortholog_TPM)
	rm(tmp_mean_zebra)
	gc()
	gc()
}


# vs Chiken
if(1){

	if(1){

		RBBH_table_chick    <- read.table("~/Documents/work/Research/Fluctuation/Data/ExpressionTables/RBBH/RBBH_gene_Oryzias_latipes_and_Gallus_gallus.txt", header=FALSE, sep="\t")
    	TPM_chick           <- read.table("/Volumes/HDD_BackUp/Fluctuation/Data/Gg_tpm.coding.txt", header=TRUE, sep="\t", row.names = 1)
		TPM_chick			<- log(TPM_chick+1, base=10)


        tmp_ortholog                        <- data.frame(RBBH_table_chick[,1])
        rownames(tmp_ortholog)              <- tmp_ortholog[,1]
        medaka_ortholog_TPM                 <- merge(tmp_inbred, tmp_ortholog, by.x = "row.names", by.y = "row.names")
        rownames(medaka_ortholog_TPM)       <- medaka_ortholog_TPM[,1]
        medaka_ortholog_TPM                 <- medaka_ortholog_TPM[,2:(ncol(medaka_ortholog_TPM)-1)]    

        TPM_chick2                          <- TPM_chick[,grep("HH16", colnames(TPM_chick))]
        tmp_ortholog                        <- data.frame(RBBH_table_chick)
        rownames(tmp_ortholog)              <- tmp_ortholog[,2]
        chick_ortholog_TPM                  <- merge(TPM_chick2, tmp_ortholog, by.x = "row.names", by.y = "row.names")
        rownames(chick_ortholog_TPM)        <- chick_ortholog_TPM[,4]
        chick_ortholog_TPM                  <- chick_ortholog_TPM[,2:3] 

        combined_RBBH                       <- merge(medaka_ortholog_TPM, chick_ortholog_TPM, by.x = "row.names", by.y = "row.names")
        rownames(combined_RBBH)             <- combined_RBBH[,1]
        combined_RBBH                       <- combined_RBBH[,-1]
        combined_RBBH                       <- na.omit(combined_RBBH)



		GXP_change_nocorrection			<- apply(combined_RBBH, 1, get_GXPchange_log)
		GXP_change_nocorrection			<- data.frame(GXP_change_nocorrection)
		tmp_gene						<- intersect(names(tmp_flu_corrected), rownames(GXP_change_nocorrection))


		tmp_explevel					<- apply(tmp_inbred[tmp_gene,], 1, mean)
		tmp_explevel					<- tmp_explevel[tmp_gene]
		tmp_explevel					<- tmp_explevel[order(tmp_explevel, decreasing=TRUE)]
		GXP_change_nocorrection			<- GXP_change_nocorrection[names(tmp_explevel),]
		names(GXP_change_nocorrection)	<- names(tmp_explevel)

		RM_GXP_change			<- numeric(length(GXP_change_nocorrection))
		names(RM_GXP_change)	<- names(GXP_change_nocorrection)
		for(k in 1:length(GXP_change_nocorrection)){
		
			if(k == 1 || k == length(GXP_change_nocorrection)){
				RM_GXP_change[k]	<- GXP_change_nocorrection[k]
			}else{
		
				if(k <= 250){
					RM_GXP_change[k]	<- median(GXP_change_nocorrection[1:(2*k-1)], na.rm=TRUE)
				}else{
					if(k <= (length(GXP_change_nocorrection) - 250)){
						RM_GXP_change[k]	<- median(GXP_change_nocorrection[(k-250):(k+250)], na.rm=TRUE)
					}else{
						RM_GXP_change[k]	<- median(GXP_change_nocorrection[(2*k-length(GXP_change_nocorrection)):length(GXP_change_nocorrection)], na.rm=TRUE)
					}						
				}					
			}	
		}

		GXP_change_corrected	<- GXP_change_nocorrection - RM_GXP_change
		GXP_change_corrected	<- data.frame(GXP_change_corrected)
    }

    if(1){
        merged              <- merge(tmp_flu_corrected, GXP_change_corrected, by.x = "row.names", by.y = "row.names") 
        rownames(merged)    <- merged[,1]
        merged              <- merged[,-1]

        outputpath  <- paste(outputdir_tmp2, "RMCorrection_fluctuation-macroevo_chick.pdf", sep="")
        pdf(file = outputpath,  width = 7,      height = 7)
		par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i")
        #par(xpd = TRUE, mar = c(8, 6, 2, 2), mgp = c(5, 1, 0), pty="s")

		color  = densCols(merged[,1], merged[,2], colramp = colorRampPalette(c("#01579B","#0277BD","#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
        plot(   merged[,1], merged[,2],
                xlab        ="Fluctuation (RM correction)", #Fluctuation
                ylab        = "Interspecies diversity (RM correction)", #Macroevolutionary change
                #log            = "xy",
                #xaxt   = "n",
                #yaxt   = "n",
				xlim		= c(-0.075, 0.25), 
				ylim		= c(-1, 2.5), 
                pch         = 16,
                col         = color, #"dodgerblue",
                cex         = 0.3,
                cex.lab     = 1,
                main        = paste("After gene selection, Log, TPM=>",TPMCutOff, ", chick", sep=""),
                las         = 1
            )
        
        mtext(paste(nrow(na.omit(merged)), " genes", sep=""), side = 1, line = 6, adj = 1)
        cortest_list    <- cor.test(merged[,1], merged[,2], method = "spearman")
        mtext(paste("fluctuation-macroevo(RM correction for both) rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
    
        dev.off()                                                                   
    }


    rm(RBBH_table_chick)
    rm(TPM_chick)
    rm(TPM_chick2)
    rm(chick_ortholog_TPM)
    rm(tmp_mean_chick)
    gc()
    gc()
}


# vs Mouse
if(1){

	if(1){
    	RBBH_table_mouse    <- read.table("~/Documents/work/Research/Fluctuation/Data/ExpressionTables/RBBH/RBBH_gene_Oryzias_latipes_and_Mus_musculus.txt", header=FALSE, sep="\t")
    	TPM_mouse           <- read.table("/Volumes/HDD_BackUp/Fluctuation/Data/Mm_tpm.coding.txt", header=TRUE, sep="\t", row.names = 1)
		TPM_mouse			<- log(TPM_mouse+1, base=10)


        tmp_ortholog                        <- data.frame(RBBH_table_mouse[,1])
        rownames(tmp_ortholog)              <- tmp_ortholog[,1]
        medaka_ortholog_TPM                 <- merge(tmp_inbred, tmp_ortholog, by.x = "row.names", by.y = "row.names")
        rownames(medaka_ortholog_TPM)       <- medaka_ortholog_TPM[,1]
        medaka_ortholog_TPM                 <- medaka_ortholog_TPM[,2:(ncol(medaka_ortholog_TPM)-1)]    

        TPM_mouse2                          <- TPM_mouse[,grep("E9.0_", colnames(TPM_mouse))]
        tmp_ortholog                        <- data.frame(RBBH_table_mouse)
        rownames(tmp_ortholog)              <- tmp_ortholog[,2]
        mouse_ortholog_TPM                  <- merge(TPM_mouse2, tmp_ortholog, by.x = "row.names", by.y = "row.names")
        rownames(mouse_ortholog_TPM)        <- mouse_ortholog_TPM[,4]
        mouse_ortholog_TPM                  <- mouse_ortholog_TPM[,2:3] 

        combined_RBBH                       <- merge(medaka_ortholog_TPM, mouse_ortholog_TPM, by.x = "row.names", by.y = "row.names")
        rownames(combined_RBBH)             <- combined_RBBH[,1]
        combined_RBBH                       <- combined_RBBH[,-1]
        combined_RBBH                       <- na.omit(combined_RBBH)


		GXP_change_nocorrection			<- apply(combined_RBBH, 1, get_GXPchange_log)
		GXP_change_nocorrection			<- data.frame(GXP_change_nocorrection)
		tmp_gene						<- intersect(names(tmp_flu_corrected), rownames(GXP_change_nocorrection))


		tmp_explevel					<- apply(tmp_inbred[tmp_gene,], 1, mean)
		tmp_explevel					<- tmp_explevel[tmp_gene]
		tmp_explevel					<- tmp_explevel[order(tmp_explevel, decreasing=TRUE)]
		GXP_change_nocorrection			<- GXP_change_nocorrection[names(tmp_explevel),]
		names(GXP_change_nocorrection)	<- names(tmp_explevel)

		RM_GXP_change			<- numeric(length(GXP_change_nocorrection))
		names(RM_GXP_change)	<- names(GXP_change_nocorrection)
		for(k in 1:length(GXP_change_nocorrection)){
		
			if(k == 1 || k == length(GXP_change_nocorrection)){
				RM_GXP_change[k]	<- GXP_change_nocorrection[k]
			}else{
			
				if(k <= 250){
					RM_GXP_change[k]	<- median(GXP_change_nocorrection[1:(2*k-1)], na.rm=TRUE)
				}else{
					if(k <= (length(GXP_change_nocorrection) - 250)){
						RM_GXP_change[k]	<- median(GXP_change_nocorrection[(k-250):(k+250)], na.rm=TRUE)
					}else{
						RM_GXP_change[k]	<- median(GXP_change_nocorrection[(2*k-length(GXP_change_nocorrection)):length(GXP_change_nocorrection)], na.rm=TRUE)
					}						
				}					
			}	
		}

		GXP_change_corrected	<- GXP_change_nocorrection - RM_GXP_change
		GXP_change_corrected	<- data.frame(GXP_change_corrected)
	}

  
    if(1){
        merged              <- merge(tmp_flu_corrected, GXP_change_corrected, by.x = "row.names", by.y = "row.names") 
        rownames(merged)    <- merged[,1]
        merged              <- merged[,-1]

        outputpath  <- paste(outputdir_tmp2, "RMCorrection_fluctuation-macroevo_mouse.pdf", sep="")
        pdf(file = outputpath,  width = 7,      height = 7)
		par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i")
        #par(xpd = TRUE, mar = c(8, 6, 2, 2), mgp = c(5, 1, 0), pty="s")

		color  = densCols(merged[,1], merged[,2], colramp = colorRampPalette(c("#01579B","#0277BD","#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
        plot(   merged[,1], merged[,2],
                xlab        ="Fluctuation (RM correction)", #Fluctuation
                ylab        = "Interspecies diversity (RM correction)", #Macroevolutionary change
                #log            = "xy",
                #xaxt   = "n",
                #yaxt   = "n",
				xlim		= c(-0.1, 0.25), 
				ylim		= c(-0.5, 2.5),  
                pch         = 16,
                col         = color, #"dodgerblue",
                cex         = 0.3,
                cex.lab     = 1,
                main        = paste("After gene selection, Log, TPM=>",TPMCutOff, ", mouse", sep=""),
                las         = 1
            )
        
        mtext(paste(nrow(na.omit(merged)), " genes", sep=""), side = 1, line = 6, adj = 1)
        cortest_list    <- cor.test(merged[,1], merged[,2], method = "spearman")
        mtext(paste("fluctuation-macroevo(RM correction for both) rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)
    
        dev.off()                                                                   
    }


    rm(RBBH_table_mouse)
    rm(TPM_mouse)
    rm(TPM_mouse2)
    rm(mouse_ortholog_TPM)
    rm(tmp_mean_mouse)
    gc()
    gc()
}

rm(GXP_change_nocorrection_conventional)
rm(GXP_change_corrected_conventional)
rm(GXP_change_nocorrection_new)
rm(GXP_change_corrected_new)
rm(merged_conventional)
rm(merged_new)
gc()
gc()


##############################################################################################################
####　Developmental stability and intraspecies diersity in Phenotype with Developmental genes			
##############################################################################################################

if(0){
	outputdir_tmp2	<- paste(outputdir_tmp, "DevelopmentalGenes/", sep="")
	dir.create(outputdir_tmp2)

	## Obtain Developmental genes("GO:0032502" developmental process, and descendants)
	if(0){
		# "GO:0032502" developmental process & descendants
		x 					<- as.list(GOBPOFFSPRING)
		developmentalgenes	<- x[["GO:0032502"]]
		write.table(devgeneGOterm, "~/Documents/work/Research/Fluctuation/Data/ExpressionTables/GO/developmentalgenes.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)


		# Get GO table
		setup_info		<- useMart(host = "jan2019.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset="olatipes_gene_ensembl")

		gene_ids		<- rownames(TPM_all)		
		ID_and_GO		<- getBM(	attributes	= c("ensembl_gene_id", "go_id"),
									filters		= "ensembl_gene_id", 
									values		= gene_ids,
									mart		= setup_info
								)

		OutputPath_tmp	<- paste(outputdir_tmp, "GOtables.txt", sep="")
		write.table(ID_and_GO,OutputPath_tmp, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)									

		tmp_devgenes	<- ID_and_GO[ID_and_GO[,2] %in% developmentalgenes,1]
		tmp_devgenes	<- unique(unlist(tmp_devgenes))

		TPM_devgenes	<- TPM_all[tmp_devgenes,]
		OutputPath_tmp	<- paste(outputdir_tmp, "DevelopmentalGenes_GO_20M.txt", sep="")
		write.table(TPM_devgenes, file = OutputPath_tmp, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)	

		TPM_devgenes	<- TPM_wild[tmp_devgenes,]
		OutputPath_tmp	<- paste(outputdir_tmp, "DevelopmentalGenes_wild_GO_20M.txt", sep="")
		write.table(TPM_devgenes, file = OutputPath_tmp, sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)	


		ID_and_GOslim	<- getBM(	attributes	= c("ensembl_gene_id", "goslim_goa_accession"),
									filters		= "ensembl_gene_id", 
									values		= gene_ids,
									mart		= setup_info
								)

		OutputPath_tmp	<- paste(outputdir_tmp,"GOSlimtables.txt", sep="")
		write.table(ID_and_GOslim, file = OutputPath_tmp, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)									

		tmp_devgenes	<- ID_and_GOslim[ID_and_GOslim[,2] %in% developmentalgenes,1]
		tmp_devgenes	<- unique(unlist(tmp_devgenes))

		TPM_devgenes	<- TPM_all[tmp_devgenes,]
		OutputPath_tmp	<- paste(outputdir_tmp,"DevelopmentalGenes_GOSlim_20M.txt", sep="")
		write.table(TPM_devgenes, file = OutputPath_tmp, sep = "\t", quote = FALSE, append=FALSE, col.names = TRUE, row.names = TRUE)	

		TPM_devgenes	<- TPM_wild[tmp_devgenes,]
		OutputPath_tmp	<- paste(outputdir_tmp,"DevelopmentalGenes_wild_GOSlim_20M.txt", sep="")
		write.table(TPM_devgenes, file = OutputPath_tmp, sep = "\t", quote = FALSE, append=FALSE, col.names = TRUE, row.names = TRUE)	
	}
	TPM_dev			<- read.table("~/Documents/work/Research/Fluctuation/Data/ExpressionTables/GOtable_developmentalgenes/DevelopmentalGenes_GO_20M.txt", header=TRUE, sep="\t", row.names=1)

	# Developmental stability
	if(1){	
		TPM_early_dev	<- TPM_early[intersect(rownames(na.omit(TPM_early)),rownames(TPM_dev)),]
		TPM_mid_dev		<- TPM_mid[intersect(rownames(na.omit(TPM_mid)),rownames(TPM_dev)),]
		TPM_late1_dev	<- TPM_late1[intersect(rownames(na.omit(TPM_late1)),rownames(TPM_dev)),]
		TPM_late2_dev	<- TPM_late2[intersect(rownames(na.omit(TPM_late2)),rownames(TPM_dev)),]	

		dist_early			<- Make_dist_array(TPM_early_dev)
		dist_mid			<- Make_dist_array(TPM_mid_dev)
		dist_late1			<- Make_dist_array(TPM_late1_dev)
		dist_late2			<- Make_dist_array(TPM_late2_dev)

		outputpath	<- paste(outputdir_tmp2, "Box_inbred_dev.pdf", sep="")
		pdf(file = outputpath,	width = 7,		height = 7)
   		par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1, dist_late2,
								horizontal	= TRUE,
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
									ylim	= c(0,0.08),
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					cex.axis	= 1,
									cex 	= 0.3,
									pch 	= 16,
									main	= paste("After gene selection, Log, inbred, dev genes with TPM=>",TPMCutOff, sep=""),
									las		= 1
								)

		par(new=T)

		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xlim	= c(0,0.08),
					xaxt		= "n",
					yaxt		= "n",
					ylab		= "",
					xlab		= "Var of TPM difference",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5)),   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
					cex.axis	= 1,
					cex.lab		= 1
					)

		axis(	side		= 4, 		
				at			= 1:4,  	
    			labels		= c(paste( nrow(TPM_early_dev), " genes",sep=""),
    							paste( nrow(TPM_mid_dev), " genes",sep=""),
    							paste( nrow(TPM_late1_dev), " genes",sep=""),
    							paste( nrow(TPM_late2_dev), " genes",sep="")
    								),
    			tck			= 0,			
    			las			= 1,  			
   		 		mgp			= c(3.5,0.7,0), 
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)

		# Statistical test Kruskal Wallis)
		KW_result <- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)

		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")
		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()
	}
	
	
	# Intraspecies diversity
	if(1){
		dist_early			<- Make_dist_array_wild(TPM_Kasasa_early[intersect(rownames(na.omit(TPM_early)),rownames(TPM_dev)),], TPM_Oura_early[intersect(rownames(na.omit(TPM_early)),rownames(TPM_dev)),])
		dist_mid			<- Make_dist_array_wild(TPM_Kasasa_mid[intersect(rownames(na.omit(TPM_mid)),rownames(TPM_dev)),], TPM_Oura_mid[intersect(rownames(na.omit(TPM_mid)),rownames(TPM_dev)),])
		dist_late1			<- Make_dist_array_wild(TPM_Kasasa_late1[intersect(rownames(na.omit(TPM_late1)),rownames(TPM_dev)),], TPM_Oura_late1[intersect(rownames(na.omit(TPM_late1)),rownames(TPM_dev)),])
		dist_late2			<- Make_dist_array_wild(TPM_Kasasa_late2[intersect(rownames(na.omit(TPM_late2)),rownames(TPM_dev)),], TPM_Oura_late2[intersect(rownames(na.omit(TPM_late2)),rownames(TPM_dev)),])

		outputpath	<- paste(outputdir_tmp2, "Box_Spearman_wild_dev.pdf", sep="")
		pdf(file = outputpath,	width = 7,		height = 7)
    	par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1, dist_late2,
									horizontal	= TRUE,
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									ylim	= c(0,0.08),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					cex.axis	= 1,
									cex 	= 0.3,
									pch 	= 16,
									main 	= paste("After gene selection, Log, wild, dev genes with TPM=>",TPMCutOff, sep=""),
									las		= 1
								)

		par(new=T)

		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim	= c(0,0.08),
					ylab		= "",
					xlab		= "Var of TPM difference",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5))   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
				)

		axis(	side		= 4, 	
				at			= 1:4,  
    			labels		= c(paste( nrow(TPM_early_dev), " genes",sep=""),
    							paste( nrow(TPM_mid_dev), " genes",sep=""),
    							paste( nrow(TPM_late1_dev), " genes",sep=""),
    							paste( nrow(TPM_late2_dev), " genes",sep="")
    								),
    			tck			= 0,				
    			las			= 1,  			
   		 		mgp			= c(3.5,0.7,0), 	
   		 		cex.axis	= 1,
   		 		adj			= 0.5
			)



		# Statistical test Kruskal Wallis)
		KW_result <- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)

		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")
				
		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()

	}

	rm(box_position)
	rm(dist)
	rm(DistStage)
	rm(stage)
	rm(dist_early)
	rm(dist_mid)
	rm(dist_late1)
	rm(dist_late2)
	rm(TPM_dev)
	
	gc()
	gc()
}




##############################################################################################################
####　Developmental stability and intraspecies diersity in Phenotype with constitutive genes			
##############################################################################################################

if(0){

	outputdir_tmp2	<- paste(outputdir_tmp, "ConstitutiveGenes/", sep="")
	dir.create(outputdir_tmp2)

	# Obtain genes with TPM≥ threshold at all stages
	if(1){
		TPM_timecourse			<- read.table("~/Documents/work/Research/Fluctuation/Data/ExpressionTables/temporal_pleiotropy/Medaka_TPM_TimeCourse.txt", header=TRUE, sep="\t", row.names=1)
		TPM_timecourse_mean		<- matrix(nrow = nrow(TPM_timecourse), ncol = 18)
		for(i in 1:18){
			TPM_timecourse_mean[,i]	<- apply(TPM_timecourse[,(3*i-2):(3*i)], 1, mean)								
		}
							
		convert_0or1	<- function(x){
			if (x >= TPMCutOff){
				return(1)
			}else{
				return(0)
			}
		}
	
		TPM_timecourse_mean				<- apply(TPM_timecourse_mean, c(1,2), convert_0or1)
		TPM_timecourse_stages			<- data.frame(apply(TPM_timecourse_mean, 1, sum))
		rownames(TPM_timecourse_stages)	<- rownames(TPM_timecourse)
	
		TPM_always_includeAdult			<- subset(TPM_timecourse_stages, TPM_timecourse_stages == 18)
		OutputPath_tmp	<- paste(outputdir_tmp2,"GenesAlwaysExpressingIncludingAdults_20M.txt", sep="")
		write.table(rownames(TPM_always_includeAdult), file = OutputPath_tmp, sep = "\t", quote = FALSE, append=FALSE, col.names = F, row.names = F)	
	}
	
	# Developmental stability
	if(1){	
		dist_early			<- Make_dist_array(TPM_early[intersect(rownames(na.omit(TPM_early)),rownames(TPM_always_includeAdult)),])
		dist_mid			<- Make_dist_array(TPM_mid[intersect(rownames(na.omit(TPM_mid)),rownames(TPM_always_includeAdult)),])
		dist_late1			<- Make_dist_array(TPM_late1[intersect(rownames(na.omit(TPM_late1)),rownames(TPM_always_includeAdult)),])
		dist_late2			<- Make_dist_array(TPM_late2[intersect(rownames(na.omit(TPM_late2)),rownames(TPM_always_includeAdult)),])

		outputpath	<- paste(outputdir_tmp2, "Box_inbred_con.pdf", sep="")
		pdf(file = outputpath,	width = 7,		height = 7)
   		par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

		box_position	<- boxplot(	dist_early, dist_mid, dist_late1,dist_late2,
									horizontal	= TRUE,
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					ylim	= c(0,0.09),
									cex 	= 0.3,
									pch 	= 16,
									main	= paste("After gene selection, Log, inbred(con genes), genes with TPM=>",TPMCutOff, sep=""),
									las		= 1
								)

		par(new=T)

		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim	= c(0,0.09),
					ylab		= "",
					xlab		= "Var of TPM difference",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5)),   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
					cex.axis	= 1,
					cex.lab		= 1
					)

		axis(	side		= 4, 	
				at			= 1:4,  
    			labels		= c(paste( nrow(na.omit(TPM_early[intersect(rownames(na.omit(TPM_early)),rownames(TPM_always_includeAdult)),])), " genes",sep=""),
    							paste( nrow(na.omit(TPM_mid[intersect(rownames(na.omit(TPM_mid)),rownames(TPM_always_includeAdult)),])), " genes",sep=""),
    							paste( nrow(na.omit(TPM_late1[intersect(rownames(na.omit(TPM_late1)),rownames(TPM_always_includeAdult)),])), " genes",sep=""),
    							paste( nrow(na.omit(TPM_late2[intersect(rownames(na.omit(TPM_late2)),rownames(TPM_always_includeAdult)),])), " genes",sep="")
    								),
    			tck			= 0,		
    			las			= 1,  		
   		 		mgp			= c(3.5,0.7,0), 	
   		 		cex.axis	= 1,
   		 		adj			= 0.5
				)



		# Statistical test Kruskal Wallis)
		KW_result <- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)

		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")
		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()

		rm(dist)
		rm(DistStage)
		rm(stage)
		gc()
		gc()

	}

	
	
	# Intraspecies diversity
	if(1){
		dist_early			<- Make_dist_array_wild(TPM_Kasasa_early[intersect(rownames(na.omit(TPM_Kasasa_early)),rownames(TPM_always_includeAdult)),], TPM_Oura_early[intersect(rownames(na.omit(TPM_Oura_early)),rownames(TPM_always_includeAdult)),])
		dist_mid			<- Make_dist_array_wild(TPM_Kasasa_mid[intersect(rownames(na.omit(TPM_Kasasa_mid)),rownames(TPM_always_includeAdult)),], TPM_Oura_mid[intersect(rownames(na.omit(TPM_Oura_mid)),rownames(TPM_always_includeAdult)),])
		dist_late1			<- Make_dist_array_wild(TPM_Kasasa_late1[intersect(rownames(na.omit(TPM_Kasasa_late1)),rownames(TPM_always_includeAdult)),], TPM_Oura_late1[intersect(rownames(na.omit(TPM_Oura_late1)),rownames(TPM_always_includeAdult)),])
		dist_late2			<- Make_dist_array_wild(TPM_Kasasa_late2[intersect(rownames(na.omit(TPM_Kasasa_late2)),rownames(TPM_always_includeAdult)),], TPM_Oura_late2[intersect(rownames(na.omit(TPM_Oura_late2)),rownames(TPM_always_includeAdult)),])

		outputpath	<- paste(outputdir_tmp2, "Box_wild_con.pdf", sep="")
		pdf(file = outputpath,	width = 7,		height = 7)
    	par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s",xaxs="i", yaxs="i")
		
		box_position	<- boxplot(	dist_early, dist_mid, dist_late1, dist_late2,
									horizontal	= TRUE,
									col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
									border	= "grey35",
					       			names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 					ylim	= c(0,0.09),
				 					cex.axis	= 1,
									cex 	= 0.3,
									pch 	= 16,
									main 	= paste("After gene selection, Log, wild(con genes), genes with TPM=>",TPMCutOff, sep=""),
									las		= 1
								)

		par(new=T)

		dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
		stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= F,
					xaxt		= "n",
					yaxt		= "n",
					xlim	= c(0,0.09),
					ylab		= "",
					xlab		= "Var of TPM difference",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
									bg=rgb(151/255,255/255,255/255, alpha=0.5),
									bg=rgb(144/255,238/255,144/255, alpha=0.5), 
									bg=rgb(255/255,193/255,37/255,  alpha=0.5)),   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
					cex.axis	= 1,
					cex.lab		= 1
				)

		axis(	side		= 4, 	
				at			= 1:4,  	
    			labels		= c(paste( nrow(TPM_Kasasa_early[intersect(rownames(na.omit(TPM_Kasasa_early)),rownames(TPM_always_includeAdult)),]), " genes",sep=""),
    							paste( nrow(TPM_Kasasa_mid[intersect(rownames(na.omit(TPM_Kasasa_early)),rownames(TPM_always_includeAdult)),]), " genes",sep=""),
    							paste( nrow(TPM_Kasasa_late1[intersect(rownames(na.omit(TPM_Kasasa_early)),rownames(TPM_always_includeAdult)),]), " genes",sep=""),
    							paste( nrow(TPM_Kasasa_late2[intersect(rownames(na.omit(TPM_Kasasa_early)),rownames(TPM_always_includeAdult)),]), " genes",sep="")
    								),
    			tck			= 0,		
    			las			= 1,  			
   		 		mgp			= c(3.5,0.7,0), 	
   		 		cex.axis	= 1,
   		 		adj			= 0.5
				)

		# Statistical test Kruskal Wallis)
		KW_result <- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
		mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)

		# Statistical test (Steel Dwass)
		stage 		<- as.factor(stage)
		SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")
		mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
		mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
		mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
		mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
		mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()

		rm(dist)
		rm(DistStage)
		rm(stage)
		gc()
		gc()

	}

	rm(TPM_timecourse)
	rm(TPM_timecourse_mean)
	rm(TPM_timecourse_stages)
	rm(TPM_always_includeAdult)	
	gc()
	gc()



}



##############################################################################################################
####　Analysis using genes in different expression ranges				
##############################################################################################################

if(1){
	outputdir_tmp2	<- paste(outputdir_tmp, "ExpressionRange/", sep="")
	dir.create(outputdir_tmp2)


	# Keep only genes with expression levels within a threshold range
		upper_lower_cut		<- function(x){
				if(sum(is.na(x))>0){
					output		<- rep(NA, length(x))
				}else{
					tmp_mean	<- mean(x)
					if(x >= lower_limit && x <= upper_limit){ 
						output 	<- unlist(x)
					 }else{
					 	 output	<- rep(NA, length(x)) 
					 }
				}
				return (output)
		}


	rangebox	<- function(percentile_low_GXP , percentile_high_GXP){


			TPM_subset				<- TPM_early			
			TPM_mean				<- apply(TPM_subset, 1, mean)			
			TPM_mean				<- na.omit(TPM_mean)
			TPM_mean				<- TPM_mean[order(TPM_mean)]			
			TPM_mean				<- TPM_mean[round(quantile(1:length(TPM_mean),(percentile_low_GXP/100))):round(quantile(1:length(TPM_mean),(percentile_high_GXP/100)))]
			TPM_percent_early		<- TPM_subset[names(TPM_mean),]

			TPM_subset				<- TPM_mid
			TPM_mean				<- apply(TPM_subset, 1, mean)			
			TPM_mean				<- na.omit(TPM_mean)
			TPM_mean				<- TPM_mean[order(TPM_mean)]			
			TPM_mean				<- TPM_mean[round(quantile(1:length(TPM_mean),(percentile_low_GXP/100))):round(quantile(1:length(TPM_mean),(percentile_high_GXP/100)))]
			TPM_percent_mid			<- TPM_subset[names(TPM_mean),]

			TPM_subset				<- TPM_late1
			TPM_mean				<- apply(TPM_subset, 1, mean)			
			TPM_mean				<- na.omit(TPM_mean)
			TPM_mean				<- TPM_mean[order(TPM_mean)]			
			TPM_mean				<- TPM_mean[round(quantile(1:length(TPM_mean),(percentile_low_GXP/100))):round(quantile(1:length(TPM_mean),(percentile_high_GXP/100)))]
			TPM_percent_late1		<- TPM_subset[names(TPM_mean),]

			TPM_subset				<- TPM_late2
			TPM_mean				<- apply(TPM_subset, 1, mean)			
			TPM_mean				<- na.omit(TPM_mean)
			TPM_mean				<- TPM_mean[order(TPM_mean)]			
			TPM_mean				<- TPM_mean[round(quantile(1:length(TPM_mean),(percentile_low_GXP/100))):round(quantile(1:length(TPM_mean),(percentile_high_GXP/100)))]
			TPM_percent_late2		<- TPM_subset[names(TPM_mean),]

			tmp_merge				<- merge(TPM_percent_late2, TPM_percent_early, by.x="row.names", by.y="row.names", all=T)		
			rownames(tmp_merge)		<- tmp_merge[,1]
			tmp_merge				<- tmp_merge[,-1]
			tmp_merge2				<- merge(tmp_merge, TPM_percent_mid, by.x="row.names", by.y="row.names", all=T)		
			rownames(tmp_merge2)	<- tmp_merge2[,1]
			tmp_merge2				<- tmp_merge2[,-1]
			tmp_merge3				<- merge(tmp_merge2, TPM_percent_late1, by.x="row.names", by.y="row.names", all=T)		
			rownames(tmp_merge3)	<- tmp_merge3[,1]
			tmp_merge3				<- tmp_merge3[,-1]
			tmp_merge				<- data.frame(gene_names[rownames(tmp_merge3)], tmp_merge3)
			colnames(tmp_merge)		<- name_inbred
			output_path_tmp			<- paste(outputdir_tmp2, "TPM_inbred-GeneNUMpercent_", percentile_low_GXP, "-", percentile_high_GXP, "_20M.txt" , sep="")
			write.table(tmp_merge, output_path_tmp, sep="\t", quote=FALSE, col.names= TRUE, row.names= TRUE, append=FALSE)


			dist_early	<- Make_dist_array(na.omit(TPM_percent_early))
			dist_mid	<- Make_dist_array(na.omit(TPM_percent_mid))
			dist_late1	<- Make_dist_array(na.omit(TPM_percent_late1))
			dist_late2	<- Make_dist_array(na.omit(TPM_percent_late2))

			outputpath	<- paste(outputdir_tmp2, "Box_inbred_", percentile_low_GXP, "-", percentile_high_GXP, ".pdf", sep="")
			pdf(file = outputpath,	width = 7,		height = 7)
   			par(xpd = TRUE, mar = c(13, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", xaxs="i", yaxs="i")

			box_position	<- boxplot(	dist_early, dist_mid, dist_late1,dist_late2,
										horizontal	= TRUE,
										col		= c("dodgerblue4","darkturquoise","darkseagreen1","lightGoldenrod1"),
										border	= "grey35",
					       				names	= c(	"Early\n(st.15)","Mid   \n(st.23.5)","Mid  \n(st.28)","Late  \n(Hatch)"),
				 						ylim	= c(0, 0.1),
				 						cex.axis	= 1,
										cex 	= 0.3,
										pch 	= 16,
										main	= paste("Log, inbred, genes with TPM=>",TPMCutOff, " ",percentile_low_GXP, "-", percentile_high_GXP,sep=""),
										las		= 1
									)

			axis(	side		= 4, 				#  右
					at			= 1:4,  		# 座標
    				labels		= c(paste( nrow(TPM_percent_early), " genes",sep=""),
    								paste( nrow(TPM_percent_mid), " genes",sep=""),
    								paste( nrow(TPM_percent_late1), " genes",sep=""),
    								paste( nrow(TPM_percent_late2), " genes",sep="")
    								),
    				tck			= 0,				#　目盛りの長さ、内側に食い込む 
    				las			= 1,  				# ラベルスタイル（las）は1（常に水平）
   		 			mgp			= c(3.5,0.7,0), 		# ラベル位置
   		 			cex.axis	= 1,
   		 			adj			= 0.5
				)

			par(new=T)
			
			dist				<- c(dist_early,dist_mid,dist_late1,dist_late2)
			stage				<- c(rep("Early",length(dist_early)),rep("Mid",length(dist_mid)),rep("Late1",length(dist_late1)),rep("Late2",length(dist_late2)))
			DistStage			<- data.frame(dist,stage)
			DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("Early","Mid","Late1","Late2"))

			stripchart(	DistStage$dist ~ DistStage$stage,
						method 		= "jitter", 
						vertical 	= F,
						xaxt		= "n",
						yaxt		= "n",
						xlim		= c(0, 0.1),
						ylab		= "",
						xlab		= "VAR of TPM difference",
						cex			= 1.5,
						pch 		= 16, 
						col 		= c(bg=rgb(0,154/255,205/255, alpha=0.5),
										bg=rgb(151/255,255/255,255/255, alpha=0.5),
										bg=rgb(144/255,238/255,144/255, alpha=0.5), 
										bg=rgb(255/255,193/255,37/255,  alpha=0.5)),   #deepskyblue3 ,DarkSlateGray1, lightgreen, lightGoldenrod1
					)

		# Statistical test Kruskal Wallis)
			KW_result	<- kruskal.test(x=list(dist_early,dist_mid,dist_late1,dist_late2))
			mtext(paste("P = ", signif(KW_result$p.value, 3), " (Kruskal Wallis)", sep=""), side = 1, line = 5, adj = 1, cex = 0.8)


		# Statistical test (Steel Dwass)
			stage 		<- as.factor(stage)
			SD_pvalues	<- pSDCFlig(dist, stage, method = "Asymptotic")	
			mtext(paste("Early-Mid(Steel-Dwass)", signif(SD_pvalues$p.val[3], digits = 3), sep = "\t"), side = 1, line = 6, adj = 1, cex = 0.8)
			mtext(paste("Early-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[1], digits = 3), sep = "\t"), side = 1, line = 7, adj = 1, cex = 0.8)
			mtext(paste("Early-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[2], digits = 3), sep = "\t"), side = 1, line = 8, adj = 1, cex = 0.8)
			mtext(paste("Mid-Late1(Steel-Dwass)", signif(SD_pvalues$p.val[5], digits = 3), sep = "\t"), side = 1, line = 9, adj = 1, cex = 0.8)
			mtext(paste("Mid-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[6], digits = 3), sep = "\t"), side = 1, line = 10, adj = 1, cex = 0.8)
			mtext(paste("Late1-Late2(Steel-Dwass)", signif(SD_pvalues$p.val[4], digits = 3), sep = "\t"), side = 1, line = 11, adj = 1, cex = 0.8)

		dev.off()

	}

	rangebox(0,20)
	rangebox(20,40)
	rangebox(40,60)
	rangebox(60,80)
	rangebox(80,100)		



	rm(dist)
	rm(DistStage)
	rm(stage)
	gc()
	gc()
	
}




##############################################################################################################
####　GO analysis		
##############################################################################################################

# GOSlim enrichment for genes with top/bottom10% developmental stability
if(1){

	outputdir_tmp2	<- paste(outputdir_tmp, "GOSlimEnrichment_Fisher-FDR_BPonly/", sep="")
	dir.create(outputdir_tmp2)

	## Sort the genes in order of developmental stability and get gene IDs.
	fluctuation_early	<- na.omit(flu_RMcorrection_early)
	fluctuation_mid		<- na.omit(flu_RMcorrection_mid)
	fluctuation_late1	<- na.omit(flu_RMcorrection_late1)
	fluctuation_late2	<- na.omit(flu_RMcorrection_late2)
	
	tmp_nrow	 		<- max(length(fluctuation_early), length(fluctuation_mid), length(fluctuation_late1), length(fluctuation_late2))

	geneID_early		<- names(fluctuation_early[order(fluctuation_early, decreasing = TRUE)])
	geneID_early		<- c(geneID_early, rep(NA, (tmp_nrow - length(geneID_early))))
	
	geneID_mid			<- names(fluctuation_mid[order(fluctuation_mid, decreasing = TRUE)])
	geneID_mid			<- c(geneID_mid, rep(NA, (tmp_nrow - length(geneID_mid))))
	
	geneID_late1		<- names(fluctuation_late1[order(fluctuation_late1, decreasing = TRUE)])
	geneID_late1		<- c(geneID_late1, rep(NA, (tmp_nrow - length(geneID_late1))))
		
	geneID_late2		<- names(fluctuation_late2[order(fluctuation_late2, decreasing = TRUE)])
	geneID_late2		<- c(geneID_late2, rep(NA, (tmp_nrow - length(geneID_late2))))

	output_path				<- paste(outputdir_tmp2, "FluctuationTable_geneID_fromLARGE.txt", sep = "")
	write.table(cbind(geneID_early,geneID_mid,geneID_late1,geneID_late2), file = output_path, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE )


	## Make annotation data ... mode: dataframe; [GO_ID, Evidence Code, Gene ID]
	setup_info		<- useMart(host = "jan2019.archive.ensembl.org", biomart = "ENSEMBL_MART_ENSEMBL", dataset="olatipes_gene_ensembl")

	early_top		<- data.frame()
	mid_top			<- data.frame()
	late1_top		<- data.frame()
	late2_top		<- data.frame()
	early_top_dev	<- data.frame()
	mid_top_dev		<- data.frame()
	late1_top_dev	<- data.frame()
	late2_top_dev	<- data.frame()	

	early_bottom	<- data.frame()
	mid_bottom		<- data.frame()
	late1_bottom	<- data.frame()
	late2_bottom	<- data.frame()
	early_bottom_dev<- data.frame()
	mid_bottom_dev	<- data.frame()
	late1_bottom_dev<- data.frame()
	late2_bottom_dev<- data.frame()


	# GO Slim P<=0.01
	stages			<- c("early", "mid", "late1", "late2")
	for(i in 1:4){
      		stage               <- stages[i]

    	    geneids             <- na.omit(get(paste("geneID_", stage, sep="")))
			result              <- getBM(   attributes  = c("go_id", "goslim_goa_accession", "ensembl_gene_id"),
                                       		filters     = "ensembl_gene_id", 
                                        	values      = geneids,
                                       		mart            = setup_info
                                    	)
       		colnames(result)    <- c("GO_id","GOSlim_id", "gene_id")        
        	annotation_data     <- as.data.frame(result)

			# Extract descentant GOSlim terms of GO:0008150 biological process
           	x                   <- as.list(GOBPOFFSPRING)
        	selectGO            <- x[["GO:0008150"]]
        	selectGO            <- c(selectGO, "GO:0008150")
			annotation_data		<- annotation_data[annotation_data[,1] %in%  selectGO,] 
			annotation_data		<- annotation_data[annotation_data[,2] %in%  selectGO,]
			annotation_data		<- annotation_data[,c(2,3)]

    

			## Compare the 10% of genes with low developmental stability (large variation) vs. the remaining  genes
			if(1){
            	percentage_10           <- round(0.1*length(geneids))
            	top10per_gene           <- geneids[1:percentage_10]
				remain_gene 			<- geneids[(percentage_10+1):length(geneids)]

				GOcount_top10per_gene		<- table(annotation_data[(annotation_data[,2] %in% top10per_gene),1])
				Totalcount_top10per_gene	<- sum(GOcount_top10per_gene) 

				GOcount_remain_gene		<- table(annotation_data[(annotation_data[,2] %in% remain_gene),1])
				Totalcount_remain_gene	<- sum(GOcount_remain_gene)

				GO_shared				<- intersect(names(GOcount_top10per_gene), names(GOcount_remain_gene))


				GOenrich_function		<- function(tmp_row){
					GOcount_top10per_eachgene	<- tmp_row[1]
					GOcount_remain_eachgene		<- tmp_row[2]

					tmp_matrix 		<- matrix(nrow=2, ncol=2)
					tmp_matrix[1,]	<- c(GOcount_top10per_eachgene, GOcount_remain_eachgene)
					tmp_matrix[2,]	<- c(Totalcount_top10per_gene, Totalcount_remain_gene)

					RelativeRatio	<- (tmp_matrix[1,1]/tmp_matrix[2,1])/(tmp_matrix[1,2]/tmp_matrix[2,2])
					Pval 			<- fisher.test(tmp_matrix, alternative="greater")$p.value

					return(c(RelativeRatio, Pval))
				}

				GOcountTable			<- cbind(GOcount_top10per_gene[GO_shared], GOcount_remain_gene[GO_shared])
				GOenrichTbale			<- t(apply(GOcountTable, 1, GOenrich_function))
				GOenrichTbale			<- cbind(GOenrichTbale, rep(NA, nrow(GOenrichTbale)))

				for(k in 1:nrow(GOenrichTbale)){
					if(is.null(GOTERM[[rownames(GOenrichTbale)[k]]])){
						GOenrichTbale[k,3]	<- NA
					}else{
						GOenrichTbale[k,3]	<- try(Term(GOTERM[[rownames(GOenrichTbale)[k]]]), silent = FALSE)
						if(class(GOenrichTbale[k,3])=="try-error"){
							GOenrichTbale[k,3]	<- NA
						}
					}
				}

				GOenrichTbale			<- cbind(GOcount_top10per_gene[rownames(GOenrichTbale)], GOcount_remain_gene[rownames(GOenrichTbale)], GOenrichTbale)
				colnames(GOenrichTbale)	<- c(paste("Count_10per(total ", Totalcount_top10per_gene, ")" ,sep=""), paste("Count_remain(total ", Totalcount_remain_gene, ")", sep=""), "Relative Ratio", "Pvalue", "Term")

	           	output_path             <- paste(outputdir_tmp2, "GOSlimenrich,Pvalue", PvalaueCutOff, "_top10fluctuation_all_", stage, ".txt", sep = "")                                           
            	write.table(GOenrichTbale, file = output_path, sep = "\t", quote = FALSE)

				GOenrichTbale			<- subset(GOenrichTbale, as.numeric(GOenrichTbale[,3])>1)

				if(nrow(GOenrichTbale)==0){
					GOenrichTbale			<- array(rep(NA, 6))
					names(GOenrichTbale)	<- c(paste("Count_10per(total ", Totalcount_top10per_gene, ")" ,sep=""), paste("Count_remain(total ", Totalcount_remain_gene, ")", sep=""), "Relative Ratio", "Pvalue", "Pvalue_FDRCorrection", "Term")
					GOenrichTbale			<- t(as.matrix(GOenrichTbale))
            		output_path             <- paste(outputdir_tmp2, "GOSlimenrich,Pvalue", PvalaueCutOff, "_top10fluctuation_selected_", stage, ".txt", sep = "")                                           
            		write.table(GOenrichTbale, file = output_path, sep = "\t", quote = FALSE)

				}else{
					if(nrow(GOenrichTbale)==1){
						tmp_GOID				<- rownames(GOenrichTbale)
						GOenrichTbale			<- t(as.matrix(c(GOenrichTbale[,1:4], p.adjust(GOenrichTbale[,4], "BH"), GOenrichTbale[,5])))
						colnames(GOenrichTbale)	<- c(colnames(GOenrichTbale)[1:4], "Pvalue_FDRCorrection", "Term")
					}else{
						GOenrichTbale			<- cbind(GOenrichTbale[,1:4], p.adjust(GOenrichTbale[,4], "BH"), GOenrichTbale[,5])
						colnames(GOenrichTbale)	<- c(colnames(GOenrichTbale)[1:4], "Pvalue_FDRCorrection", colnames(GOenrichTbale)[6])
					}
            		output_path             <- paste(outputdir_tmp2, "GOSlimenrich,Pvalue", PvalaueCutOff, "_top10fluctuation_selected_", stage, ".txt", sep = "")                                        
            		write.table(GOenrichTbale, file = output_path, sep = "\t", quote = FALSE)
				}

				tmp             		<- subset(GOenrichTbale, as.numeric(GOenrichTbale[,5]) <= PvalaueCutOff)
				if(nrow(tmp)==0){
					assign((paste(stage, "_top", sep="")), matrix(nrow=0, ncol=3))
				}else{
					if(nrow(GOenrichTbale)==1){	
						assign((paste(stage, "_top", sep="")), cbind(tmp_GOID, tmp[,3], tmp[,6]))
					}else{
						assign((paste(stage, "_top", sep="")), cbind(rownames(tmp), tmp[,3], tmp[,6]))
					}
				}

				tmp_colnames	<- c("GOID", paste("Nfolds_", stage, sep=""), "Term")
           		eval(parse(text = paste("colnames(", paste(stage, "_top", sep=""), ") <- tmp_colnames", sep="")))
			}



       		## Compare the 10% of genes with high developmental stability (small variation) vs. the remaining  genes
			if(1){
            	bottom10per_gene	<- geneids[(length(geneids) - percentage_10):length(geneids)]
				remain_gene			<- geneids[1:((length(geneids) - percentage_10)-1)]

				GOcount_bottom10per_gene	<- table(annotation_data[(annotation_data[,2] %in% bottom10per_gene),1])
				Totalcount_bottom10per_gene	<- sum(GOcount_bottom10per_gene) 

				GOcount_remain_gene		<- table(annotation_data[(annotation_data[,2] %in% remain_gene),1])
				Totalcount_remain_gene	<- sum(GOcount_remain_gene)

				GO_shared				<- intersect(names(GOcount_bottom10per_gene), names(GOcount_remain_gene))


				GOenrich_function		<- function(tmp_row){
					GOcount_bottom10per_eachgene	<- tmp_row[1]
					GOcount_remain_eachgene			<- tmp_row[2]

					tmp_matrix 		<- matrix(nrow=2, ncol=2)
					tmp_matrix[1,]	<- c(GOcount_bottom10per_eachgene, GOcount_remain_eachgene)
					tmp_matrix[2,]	<- c(Totalcount_bottom10per_gene, Totalcount_remain_gene)

					RelativeRatio	<- (tmp_matrix[1,1]/tmp_matrix[2,1])/(tmp_matrix[1,2]/tmp_matrix[2,2])
					Pval 			<- fisher.test(tmp_matrix, alternative="greater")$p.value

					return(c(RelativeRatio, Pval))
				}

				GOcountTable			<- cbind(GOcount_bottom10per_gene[GO_shared], GOcount_remain_gene[GO_shared])
				GOenrichTbale			<- t(apply(GOcountTable, 1, GOenrich_function))
				GOenrichTbale			<- cbind(GOenrichTbale, rep(NA, nrow(GOenrichTbale)))				

				for(k in 1:nrow(GOenrichTbale)){
					if(is.null(GOTERM[[rownames(GOenrichTbale)[k]]])){
						GOenrichTbale[k,3]	<- NA
					}else{
						GOenrichTbale[k,3]	<- try(Term(GOTERM[[rownames(GOenrichTbale)[k]]]), silent = FALSE)
						if(class(GOenrichTbale[k,3])=="try-error"){
							GOenrichTbale[k,3]	<- NA
						}
					}
				}

				GOenrichTbale			<- cbind(GOcount_bottom10per_gene[rownames(GOenrichTbale)], GOcount_remain_gene[rownames(GOenrichTbale)], GOenrichTbale)
				colnames(GOenrichTbale)	<- c(paste("Count_10per(total ", Totalcount_bottom10per_gene, ")" ,sep=""), paste("Count_remain(total ", Totalcount_remain_gene, ")", sep=""), "Relative Ratio", "Pvalue", "Term")

	           	output_path             <- paste(outputdir_tmp2, "GOSlimenrich,Pvalue_bottom10fluctuation_all_", stage, ".txt", sep = "")                                           
            	write.table(GOenrichTbale, file = output_path, sep = "\t", quote = FALSE)

				GOenrichTbale			<- subset(GOenrichTbale, as.numeric(GOenrichTbale[,3])>1)

				if(nrow(GOenrichTbale)==0){
					GOenrichTbale			<- array(rep(NA, 6))
					names(GOenrichTbale)	<- c(paste("Count_10per(total ", Totalcount_bottom10per_gene, ")" ,sep=""), paste("Count_remain(total ", Totalcount_remain_gene, ")", sep=""), "Relative Ratio", "Pvalue", "Pvalue_FDRCorrection", "Term")
					GOenrichTbale			<- t(as.matrix(GOenrichTbale))
            		output_path             <- paste(outputdir_tmp2, "GOSlimenrich,Pvalue_bottom10fluctuation_selected_", stage, ".txt", sep = "")                                           
            		write.table(GOenrichTbale, file = output_path, sep = "\t", quote = FALSE)

				}else{
					if(nrow(GOenrichTbale)==1){
						tmp_GOID				<- rownames(GOenrichTbale)
						GOenrichTbale			<- t(as.matrix(c(GOenrichTbale[,1:4], p.adjust(GOenrichTbale[,4], "BH"), GOenrichTbale[,5])))
						colnames(GOenrichTbale)	<- c(colnames(GOenrichTbale)[1:4], "Pvalue_FDRCorrection", "Term")
						rownames(GOenrichTbale)	<- tmp_GOID
					}else{
						GOenrichTbale			<- cbind(GOenrichTbale[,1:4], p.adjust(GOenrichTbale[,4], "BH"), GOenrichTbale[,5])
						colnames(GOenrichTbale)	<- c(colnames(GOenrichTbale)[1:4], "Pvalue_FDRCorrection", colnames(GOenrichTbale)[6])
					}
            		output_path             <- paste(outputdir_tmp2, "GOSlimenrich,Pvalue_bottom10fluctuation_selected_", stage, ".txt", sep = "")                                        
            		write.table(GOenrichTbale, file = output_path, sep = "\t", quote = FALSE)
				}

				tmp             		<- subset(GOenrichTbale, as.numeric(GOenrichTbale[,5]) <= PvalaueCutOff)
				if(nrow(tmp)==0){
					assign((paste(stage, "_bottom", sep="")), matrix(nrow=0, ncol=3))
				}else{
					if(nrow(GOenrichTbale)==1){	
						assign((paste(stage, "_bottom", sep="")), cbind(tmp_GOID, tmp[,3], tmp[,6]))
					}else{
						assign((paste(stage, "_bottom", sep="")), cbind(rownames(tmp), tmp[,3], tmp[,6]))
					}
				}

				tmp_colnames	<- c("GOID", paste("Nfolds_", stage, sep=""), "Term")
           		eval(parse(text = paste("colnames(", paste(stage, "_bottom", sep=""), ") <- tmp_colnames", sep="")))



				# GO:0032502 developmental process
					tmp				<- annotation_data[annotation_data[,2] %in% bottom10per_gene,]
					tmp				<- unique(subset(tmp, tmp[,1]=="GO:0032502"))
					output_path		<- paste(outputdir_tmp2, "GO=0032502_developmental-process_bottom10fluctuation_", stage, ".txt", sep = "")                                           
            		write.table(tmp, file = output_path, sep = "\t", quote = FALSE, row.names=F, col.names=F)


				# GO:0048856 anatomical structure development
					tmp				<- annotation_data[annotation_data[,2] %in% bottom10per_gene,]
					tmp				<- unique(subset(tmp, tmp[,1]=="GO:0048856"))
					output_path		<- paste(outputdir_tmp2, "GO=0048856_anatomical-structure-development_bottom10fluctuation_", stage, ".txt", sep = "")                                           
            		write.table(tmp, file = output_path, sep = "\t", quote = FALSE, row.names=F, col.names=F)

				# GO:0048646 anatomical structure formation involved in morphogenesis
					tmp				<- annotation_data[annotation_data[,2] %in% bottom10per_gene,]
					tmp				<- unique(subset(tmp, tmp[,1]=="GO:0048646"))
					output_path		<- paste(outputdir_tmp2, "GO=0048646_anatomical-structure-formation-involved-in-morphogenesis_bottom10fluctuation_", stage, ".txt", sep = "")                                           
            		write.table(tmp, file = output_path, sep = "\t", quote = FALSE, row.names=F, col.names=F)

				# GO:0007267 cell-cell signaling
					tmp				<- annotation_data[annotation_data[,2] %in% bottom10per_gene,]
					tmp				<- unique(subset(tmp, tmp[,1]=="GO:0007267"))
					output_path		<- paste(outputdir_tmp2, "GO=0007267_cell-cell-signaling_bottom10fluctuation_", stage, ".txt", sep = "")                                           
            		write.table(tmp, file = output_path, sep = "\t", quote = FALSE, row.names=F, col.names=F)

			}


	
	
			## Plot developmental stability and intraspecies diversity with gene names
			if(1){

				tmp_top_dev		<- gene_names[unique(go_and_id_top10_dev[,2])]
	           	output_path		<- paste(outputdir_tmp2, "DevGenes-in-top10fluctuation_", stage, ".txt", sep = "")                                           
            	write.table(tmp_top_dev, file = output_path, sep = "\t", quote = FALSE, col.names=F)
				
				for(k in 1:length(tmp_top_dev)){
					if(tmp_top_dev[k] == "-"){
						tmp_top_dev[k] <- names(tmp_top_dev)[k]
					}
				}


				tmp_bottom_dev		<- gene_names[unique(go_and_id_bottom10_dev[,2])]
	           	output_path		<- paste(outputdir_tmp2, "DevGenes-in-bottom10fluctuation_", stage, ".txt", sep = "")                                           
            	write.table(tmp_bottom_dev, file = output_path, sep = "\t", quote = FALSE, col.names=F)
				
				for(k in 1:length(tmp_bottom_dev)){
					if(tmp_bottom_dev[k] == "-"){
						tmp_bottom_dev[k] <- names(tmp_bottom_dev)[k]
					}
				}

				flu_corrected		<- get(paste("flu_RMcorrection_", stage, sep=""))
				microevo_corrected	<- get(paste("microevo_RMcorrection_new_", stage, sep=""))
				tmp_gene			<- intersect(names(flu_corrected), names(microevo_corrected))
				flu_corrected		<- flu_corrected[tmp_gene]
				microevo_corrected	<- microevo_corrected[tmp_gene]



				### 描画
				if(1){
					outputdir_tmp2	<- paste(outputdir_tmp, "fluctuation,microevo,explevel,techerror_AfterGeneSelection/", sep="")
					dir.create(outputdir_tmp2)
					outputpath		<- paste(outputdir_tmp2, "RMCorrection_fluctuation-microevo(new)_devgenes-shown_withName_", stage, ".pdf", sep="")
					pdf(file = outputpath, width = 7,	height = 7)
					par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 
																												
					plot(	flu_corrected, microevo_corrected, 
							xlim	= c(-0.2, 0.8), 
							ylim	= c(-0.2, 1.2), 
							xlab	= "Fluctuation (RM correction)",
							ylab	= "Microevolutionary diversity (RM correction, new)",
							pch 	= 20,
							cex		= 0.2,
							col		= "gray60",
							cex.lab = 1,
							main	= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
							las		=1
						)
			
					par(new=T)

					plot(	flu_corrected[names(flu_corrected) %in% names(tmp_top_dev)], microevo_corrected[names(microevo_corrected) %in% names(tmp_top_dev)], 
							xlim	= c(-0.2, 0.8), 
							ylim	= c(-0.2, 1.2),
							xaxt	= "n",
							yaxt	= "n",
							xlab	= "",
							ylab	= "",
							pch 	= 20,
							cex		= 0.5,
							col		= "deeppink3",
						)

					par(new=T)

					plot(	flu_corrected[names(flu_corrected) %in% names(tmp_bottom_dev)], microevo_corrected[names(microevo_corrected) %in% names(tmp_bottom_dev)], 
							xlim	= c(-0.2, 0.8), 
							ylim	= c(-0.2, 1.2),
							xaxt	= "n",
							yaxt	= "n",
							xlab	= "",
							ylab	= "",
							pch 	= 20,
							cex		= 0.5,
							col		= "deepskyblue",
						)

					if(length(tmp_top_dev)>0){	
						text(	flu_corrected[names(flu_corrected) %in% names(tmp_top_dev)], microevo_corrected[names(microevo_corrected) %in% names(tmp_top_dev)],
						labels=paste(tmp_top_dev[names(flu_corrected[names(flu_corrected) %in% names(tmp_top_dev)])], "\n", sep=""), col="deeppink", cex=0.1)
					}

					if(length(tmp_bottom_dev)>0){	
						text(	flu_corrected[names(flu_corrected) %in% names(tmp_bottom_dev)], microevo_corrected[names(microevo_corrected) %in% names(tmp_bottom_dev)],
						labels=paste(tmp_bottom_dev[names(flu_corrected[names(flu_corrected) %in% names(tmp_bottom_dev)])], "\n", sep=""), col="dodgerblue3", cex=0.1)
					}

					dev.off()
				}

				if(1){
					outputdir_tmp2	<- paste(outputdir_tmp, "fluctuation,microevo,explevel,techerror_AfterGeneSelection/", sep="")
					dir.create(outputdir_tmp2)
					outputpath		<- paste(outputdir_tmp2, "RMCorrection_fluctuation-microevo(new)_devgenes-shown_withoutName_", stage, ".pdf", sep="")
					pdf(file = outputpath, width = 7,	height = 7)
					par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i") # 
																												
					plot(	flu_corrected, microevo_corrected, 
							xlim	= c(-0.2, 0.8), 
							ylim	= c(-0.2, 1.2), 
							xlab	= "Fluctuation (RM correction)",
							ylab	= "Microevolutionary diversity (RM correction, new)",
							pch 	= 20,
							cex		= 0.2,
							col		= "gray60",
							cex.lab = 1,
							main	= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
							las		=1
						)
			
					par(new=T)

					plot(	flu_corrected[names(flu_corrected) %in% names(tmp_top_dev)], microevo_corrected[names(microevo_corrected) %in% names(tmp_top_dev)], 
							xlim	= c(-0.2, 0.8), 
							ylim	= c(-0.2, 1.2),
							xaxt	= "n",
							yaxt	= "n",
							xlab	= "",
							ylab	= "",
							pch 	= 20,
							cex		= 0.5,
							col		= "deeppink3",
						)

					par(new=T)

					plot(	flu_corrected[names(flu_corrected) %in% names(tmp_bottom_dev)], microevo_corrected[names(microevo_corrected) %in% names(tmp_bottom_dev)], 
							xlim	= c(-0.2, 0.8), 
							ylim	= c(-0.2, 1.2),
							xaxt	= "n",
							yaxt	= "n",
							xlab	= "",
							ylab	= "",
							pch 	= 20,
							cex		= 0.5,
							col		= "deepskyblue",
						)


					dev.off()
				}

				outputdir_tmp2	<- paste(outputdir_tmp, "GOSlimEnrichment_Fisher-FDR_BPonly/", sep="")
			}	

	}



	## Heatmap (small developmental stability)
	if(1){

		tmp		<- apply(merge(early_top, mid_top, key="GOID",all=TRUE), c(1,2), as.character)
		tmp2	<- apply(merge(tmp,late1_top, key="GOID", all=TRUE), c(1,2), as.character)
		tmp3	<- merge(tmp2, late2_top, key="GOID", all=TRUE)

		## omit GO number
		tmp4	<- tmp3[,-1]		
		tmp4	<- tmp4[!is.na(tmp4[,1]),]	

		## omit Terms
		rownames(tmp4)	<- tmp4[,1]	
		tmp4			<- tmp4[,-1]
		#tmp4			<- as.matrix(tmp4)

		## character -> numeric
		tmp5 <-  apply(tmp4, c(1,2), as.numeric)

		colnames(tmp5)	<- c("Early", "Mid(st.23.5)", "Mid(st.28)", "Late(Hatch)")
			output_path	<- paste(outputdir_tmp2, "GOSlim_heatmap_HighFluctuation_P", PvalaueCutOff, ".txt", sep = "")
			write.table(tmp5, file = output_path, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE )
		colnames(tmp5)	<- c(paste("Early(st.15)\n",nrow(early_top), "genes", sep=""), paste("Mid(st.23.5)\n", nrow(mid_top), "genes", sep=""), paste("Mid(st.28)\n", nrow(late1_top), "genes", sep=""), paste("Late(Hatch)\n", nrow(late2_top), "genes", sep=""))


		output_path		<- paste(outputdir_tmp2, "GOSlim_heatmap_HighFluctuation_P", PvalaueCutOff, ".pdf", sep="")
		pdf(file = output_path, width = 14, height = 14)
		par(mar = c(5, 3, 4, 6), oma = c(1,1,1,15), xpd = TRUE) #, oma = c(4,4,4,20), mgp = c(3.5, 1, 0)

		colors <- rev(heat.colors((length(2:ceiling(max(tmp5,na.rm=TRUE)))+1)))
	
		heatmap.2(	tmp5, 
					trace			= "none",
					density.info	= "none",
					Colv			= "FALSE",
					Rowv			= "FALSE",
					cexRow			= 0.9,
					cexCol			= 0.5,
					col				= colors,
					na.color		= "gray",
					margins			= c(5,18),
					lwid			= c(2,10), 
					lhei			= c(2,15)
				)

			dev.off()
	}


	## Heatmap  (high developmental stability)
	if(1){
		tmp		<- apply(merge(early_bottom, mid_bottom, key="GOID",all=TRUE), c(1,2), as.character)
		tmp2	<- apply(merge(tmp,late1_bottom, key="GOID", all=TRUE), c(1,2), as.character)
		tmp3	<- merge(tmp2, late2_bottom, key="GOID", all=TRUE)

		## omit GO number
		tmp4	<- tmp3[,-1]		
		tmp4	<- tmp4[!is.na(tmp4[,1]),]	

		## omit Terms
		rownames(tmp4)	<- tmp4[,1]	
		tmp4			<- tmp4[,-1]
		tmp4			<- as.matrix(tmp4)

		## character -> numeric 
		tmp5 <-  apply(tmp4, c(1,2), as.numeric)

		colnames(tmp5)	<- c("Early", "Mid(st.23.5)", "Mid(st.28)", "Late(Hatch)")
			output_path	<- paste(outputdir_tmp2, "GOSlim_heatmap_LowFluctuation_P", PvalaueCutOff, ".txt", sep = "")
			write.table(tmp5, file = output_path, quote = FALSE, sep = "\t", col.names = TRUE, row.names = TRUE )
		colnames(tmp5)	<- c(paste("Early(st.15)\n",nrow(early_bottom), "genes", sep=""), paste("Mid(st.23.5)\n", nrow(mid_bottom), "genes", sep=""), paste("Mid(st.28)\n", nrow(late1_bottom), "genes", sep=""), paste("Late(Hatch)\n", nrow(late2_bottom), "genes", sep=""))


		output_path		<- paste(outputdir_tmp2, "GOSlim_heatmap_LowFluctuation_P", PvalaueCutOff, ".pdf", sep="")
		pdf(file = output_path, width = 14, height = 14)
		par(mar = c(5, 3, 4, 6), oma = c(1,1,1,15), xpd = TRUE) #, oma = c(4,4,4,20), mgp = c(3.5, 1, 0)

		colors <- rev(colorRampPalette(c("#084594", "#6baed6", "#deebf7", "#f7fbff"))(length(2:ceiling(max(tmp5,na.rm=TRUE)))+1))
	
		heatmap.2(	tmp5, 
					trace			= "none",
					density.info	= "none",
					Colv			= "FALSE",
					Rowv			= "FALSE",
					cexRow			= 0.9,
					cexCol			= 0.5,
					col				= colors,
					na.color		= "gray",
					margins			= c(5,18),
					lwid			= c(2,10), 
					lhei			= c(2,15)
				)

			dev.off()
	}



	rm(early_top)
	rm(mid_top)
	rm(late1_top)
	rm(late2_top)
	rm(early_top_dev)
	rm(mid_top_dev)
	rm(late1_top_dev)
	rm(late2_top_dev)

	rm(early_bottom)
	rm(mid_bottom)
	rm(late1_bottom)
	rm(late2_bottom)
	rm(early_bottom_dev)
	rm(mid_bottom_dev)
	rm(late1_bottom_dev)
	rm(late2_bottom_dev)

	rm(fluctuation_early)
	rm(fluctuation_mid)
	rm(fluctuation_late1)
	rm(fluctuation_late2)

	rm(geneID_early)
	rm(geneID_mid)
	rm(geneID_late1)
	rm(geneID_late2)

	rm(tmp)
	rm(tmp2)	
	rm(tmp3)
	rm(tmp4)
	rm(tmp5)

	rm(tmp_nrow)
	
	gc()
	gc()
}




##############################################################################################################
#### Correlation with cis-regulatory features
##############################################################################################################

if(1){
	outputdir_tmp2	<- paste(outputdir_tmp, "Enhancer_Features/", sep="")
	dir.create(outputdir_tmp2)

	stages		<- c("early", "mid", "late1", "late2")
	stage_num	<- c("15", "24", "28", "40")
	
	# Extract open regions from j bp upstream to k bp downstream of the TSS. Obtain the number and length of open regions and further plot.
	get_EnhancerInfo_PlotFigures	<- function(upper, lower){
		outputdir_tmp3	<- paste(outputdir_tmp2, "upper", upper, "-", "lower", lower, "/", sep="")
		dir.create(outputdir_tmp3)

		for(i in 1:4){
			stage	<- stages[i]
			
			if(1){
				# A table with information on the start point of open regions, end point of open regions, TSS, and gene ID , in that order.
				enhancer_min_bed	<- read.table(paste("/Volumes/HDD_BackUp/Fluctuation/Data/Enhancer_data/Ol_st", stage_num[i], "_OpenRegion-ClosestTSS_minimum.bed", sep=""), header=FALSE, sep="\t")

				region_withn_upperlower	<- function(enhancer_bed_row){
					openregion_center	<- (as.numeric(enhancer_bed_row[1])+as.numeric(enhancer_bed_row[2]))/2
					if(openregion_center > (as.numeric(enhancer_bed_row[3])-upper) &&  openregion_center < (as.numeric(enhancer_bed_row[3])+lower)){
						tmp_length	<- as.numeric(enhancer_bed_row[2]) - as.numeric(enhancer_bed_row[1])
						return(c(enhancer_bed_row[1], enhancer_bed_row[2], tmp_length, enhancer_bed_row[3], enhancer_bed_row[4]))
					}else{
						return(rep(NA,5))
					}
				}

				enhancer_table	<- t(apply(enhancer_min_bed, 1, region_withn_upperlower))
				enhancer_table	<- data.frame(na.omit(enhancer_table))
				colnames(enhancer_table)		<- c("start", "end", "length", "TSS", "gene_id")
				enhancer_table	<- enhancer_table[order(enhancer_table[,5]),]

				output_path_tmp = paste(outputdir_tmp3, "extracted_OpenRegion_", stage, ".text", sep="")
				write.table(enhancer_table, file = output_path_tmp, row.names=FALSE, col.names=TRUE, quote=FALSE)
			}


			# Obtain the number and length of open regions for each gene.
			if(1){
				num_length_enhancer				<- matrix(nrow=length(table(as.character(enhancer_table[,5]))), ncol=4)
				colnames(num_length_enhancer)	<- c("gene_id", "num_enhancer", "total_length_enhancer", "total_length_enhancer_withn_upperlower")

				num_length_enhancer[,1]			<- names(table(as.character(enhancer_table[,5])))

				num_length_enhancer[,2]		<- table(as.character(enhancer_table[,5]))

				counter <- 0
				for (m in 1:nrow(num_length_enhancer)){
					# enhancers of gene i: row [counter+1]~[counter + num of enhancers]
					counter_tmp 				<- counter + as.numeric(num_length_enhancer[m,2])	
					num_length_enhancer[m,3]	<- sum(as.numeric(enhancer_table[(counter+1):counter_tmp,3]))

					get_overlaplength			<- function(row){
						row	<- as.numeric(unlist(row[1:4]))
						tmp <- length(intersect(row[1]:row[2], (row[4]-upper):(row[4]+lower)))
						return(tmp)
					}
					tmp_overlaplength			<- apply(enhancer_table[(counter+1):counter_tmp,], 1, get_overlaplength)
					num_length_enhancer[m,4]	<- sum(tmp_overlaplength)
	
					counter 					<- counter_tmp		
				}

				rownames(num_length_enhancer)		<- num_length_enhancer[,1]
				num_length_enhancer					<- num_length_enhancer[,-1]
				num_length_enhancer					<- apply(num_length_enhancer, c(1,2), as.numeric)

				output_path_tmp = paste(outputdir_tmp3, "Num-Totallength-Meanlength_OpenRegion_", stage, ".text", sep="")
				write.table(num_length_enhancer, file = output_path_tmp, row.names=TRUE, col.names=TRUE, quote=FALSE)
			}


			# Visualization (vs developmental stability)
			if(1){
				tmp_flu_corrected			<- get(paste("flu_RMcorrection_", stage, sep=""))
				tmp_gene					<- intersect(rownames(num_length_enhancer), names(tmp_flu_corrected))

				outputpath	<- paste(outputdir_tmp3, "fluctuation-NumOpenRegion_", stage, ".pdf", sep="")
				pdf(file = outputpath,	width = 7,		height = 7)
				par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i")

				color  = densCols(num_length_enhancer[tmp_gene,1], tmp_flu_corrected[tmp_gene], colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
				plot(	num_length_enhancer[tmp_gene,1], tmp_flu_corrected[tmp_gene],
						xlim	= c(0, max(num_length_enhancer[tmp_gene,1])), 
						xlab		= "Number of open regions",
						ylab		= "Fluctuation (RM correction)",
						pch			= 16,
						col			= color,
						cex			= 0.6,
						main		= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
						las			= 1
					)

				mtext(paste(length(tmp_gene), " genes", sep=""), side = 1, line = 6, adj = 1)
				cortest_list	<- cor.test(num_length_enhancer[tmp_gene,1], tmp_flu_corrected[tmp_gene], method = "spearman")
				mtext(paste("fluctuation-Number of open regions rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)

				dev.off()



				outputpath	<- paste(outputdir_tmp3, "fluctuation-TotalLengthOpenRegion_", stage, ".pdf", sep="")
				pdf(file = outputpath,	width = 7,		height = 7)
				par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i")

				color  = densCols(num_length_enhancer[tmp_gene,3], tmp_flu_corrected[tmp_gene], colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
				plot(	num_length_enhancer[tmp_gene,3], tmp_flu_corrected[tmp_gene],
						xlab		= "Total length of open regions(overlapped region)",
						ylab		= "Fluctuation (RM correction)",
						log			= "x",
						pch			= 16,
						col			= color,
						cex			= 0.6,
						main		= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
						las			= 1
					)

				mtext(paste(length(tmp_gene), " genes", sep=""), side = 1, line = 6, adj = 1)
				cortest_list	<- cor.test(num_length_enhancer[tmp_gene,3], tmp_flu_corrected[tmp_gene], method = "spearman")
				mtext(paste("fluctuation-Total length of open regions rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)

				dev.off()


			}

			# Visualization (vs intraspecies diversity)
			if(1){
				tmp_microevo_corrected_new	<- get(paste("microevo_RMcorrection_new_", stage, sep=""))
				tmp_gene					<- intersect(rownames(num_length_enhancer), names(tmp_microevo_corrected_new))

				outputpath	<- paste(outputdir_tmp3, "microevo-NumOpenRegion_", stage, ".pdf", sep="")
				pdf(file = outputpath,	width = 7,		height = 7)
				par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i")

				color  = densCols(num_length_enhancer[tmp_gene,1], tmp_microevo_corrected_new[tmp_gene], colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
				plot(	num_length_enhancer[tmp_gene,1], tmp_microevo_corrected_new[tmp_gene],
						xlim	= c(0, max(num_length_enhancer[tmp_gene,1])), 
						xlab		= "Number of open regions",
						ylab		= "Microevolutionary diversity (RM correction, new)",
						pch			= 16,
						col			= color,
						cex			= 0.6,
						main		= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
						las			= 1
					)

				mtext(paste(length(tmp_gene), " genes", sep=""), side = 1, line = 6, adj = 1)
				cortest_list	<- cor.test(num_length_enhancer[tmp_gene,1], tmp_microevo_corrected_new[tmp_gene], method = "spearman")
				mtext(paste("microevo-Number of open regions rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)

				dev.off()



				outputpath	<- paste(outputdir_tmp3, "microevo-TotalLengthOpenRegion_", stage, ".pdf", sep="")
				pdf(file = outputpath,	width = 7,		height = 7)
				par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i")

				color  = densCols(num_length_enhancer[tmp_gene,3], tmp_microevo_corrected_new[tmp_gene], colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
				plot(	num_length_enhancer[tmp_gene,3], tmp_microevo_corrected_new[tmp_gene],
						xlab		= "Total length of open regions (overlapped region)",
						ylab		= "Microevolutionary diversity (RM correction, new)",
						log			= "x",
						pch			= 16,
						col			= color,
						cex			= 0.6,
						main		= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
						las			= 1
					)

				mtext(paste(length(tmp_gene), " genes", sep=""), side = 1, line = 6, adj = 1)
				cortest_list	<- cor.test(num_length_enhancer[tmp_gene,3], tmp_microevo_corrected_new[tmp_gene], method = "spearman")
				mtext(paste("microevo-Total length of open regions rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)

				dev.off()
			}

		}
	}

	get_EnhancerInfo_PlotFigures(10000, 10000)
	get_EnhancerInfo_PlotFigures(5000, 5000)
	get_EnhancerInfo_PlotFigures(3000, 3000)
	get_EnhancerInfo_PlotFigures(100, 50)
}


##############################################################################################################
#### Correlation with SNP numbers in cis-regulatory region
##############################################################################################################

if(1){
	outputdir_tmp2	<- paste(outputdir_tmp, "SNP_Features/", sep="")
	dir.create(outputdir_tmp2)


	stages		<- c("early", "mid", "late1", "late2")
	stage_num	<- c("15", "24", "28", "40")
	
	# Extract SNPs from j bp upstream to k bp downstream of the TSS. Obtain the number and length of open regions and further plot.
	get_SNPInfo_PlotFigures	<- function(upper, lower){
		outputdir_tmp3	<- paste(outputdir_tmp2, "upper", upper, "-", "lower", lower, "/", sep="")
		dir.create(outputdir_tmp3)

		for(i in 1:4){

			stage	<- stages[i]

			## A table with information on the SNP sites, TSS, and gene ID , in that order.
			SNP_TSS_bed	<- read.table(paste("/Volumes/HDD_BackUp/Fluctuation/Data/SNP_data/SNP-ClosestTSS_minimum_st", stage_num[i], ".bed", sep=""), header=FALSE, sep="\t")

			SNP_withn_upperlower	<- function(SNP_TSS_bed_bed_row){
				SNP_site	<- as.numeric(SNP_TSS_bed_bed_row[1])
				upper		<- as.numeric(SNP_TSS_bed_bed_row[2])-upper
				lower		<- as.numeric(SNP_TSS_bed_bed_row[2])+lower
				if(SNP_site > upper &&  SNP_site < lower){
					return(SNP_TSS_bed_bed_row)
				}else{
					return(rep(NA,3))
				}
			}


			SNP_table			<- t(apply(SNP_TSS_bed, 1, SNP_withn_upperlower))
			SNP_table			<- data.frame(na.omit(SNP_table))
			colnames(SNP_table)	<- c("SNP site", "TSS", "gene_id")
			SNP_table			<- unique(SNP_table[order(SNP_table[,3]),])

			output_path_tmp = paste(outputdir_tmp3, "extracted_SNP-in-OpenRegion_", stage, ".text", sep="")
			write.table(SNP_table, file = output_path_tmp, row.names=FALSE, col.names=TRUE, quote=FALSE)


			num_SNP		<- table(as.character(SNP_table[,3]))
				
			output_path_tmp = paste(outputdir_tmp3, "Num-SNP_OpenRegion_", stage, ".text", sep="")
			write.table(data.frame(num_SNP), file = output_path_tmp, row.names=FALSE, col.names=FALSE, quote=FALSE)



			# Visualization (vs developmental stability)
			if(1){
				tmp_flu_corrected			<- get(paste("flu_RMcorrection_", stage, sep=""))
				tmp_gene					<- intersect(names(num_SNP), names(tmp_flu_corrected))

				outputpath	<- paste(outputdir_tmp3, "fluctuation-Num-SNP-inOpenRegion_", stage, ".pdf", sep="")
				pdf(file = outputpath,	width = 7,		height = 7)
				par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i")

				color  = densCols(num_SNP[tmp_gene], tmp_flu_corrected[tmp_gene], colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
				plot(	as.numeric(num_SNP[tmp_gene]), as.numeric(tmp_flu_corrected[tmp_gene]),
						xlim	= c(0, (max(num_SNP[tmp_gene]+1))), 
						xlab		= "Number of substitutions in open regions",
						ylab		= "Fluctuation (RM correction)",
						pch			= 16,
						col			= color,
						cex			= 0.6,
						main		= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
						las			= 1
					)

				mtext(paste(length(tmp_gene), " genes", sep=""), side = 1, line = 6, adj = 1)
				cortest_list	<- cor.test(num_SNP[tmp_gene], tmp_flu_corrected[tmp_gene], method = "spearman")
				mtext(paste("fluctuation-Number of SNP rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)

				dev.off()

			}

			# Visualization (vs intraspecies diversity)
			if(1){
				tmp_microevo_corrected_new	<- get(paste("microevo_RMcorrection_new_", stage, sep=""))
				tmp_gene					<- intersect(names(num_SNP), names(tmp_microevo_corrected_new))

				outputpath	<- paste(outputdir_tmp3, "microevo-Num-SNP-inOpenRegion_", stage, ".pdf", sep="")
				pdf(file = outputpath,	width = 7,		height = 7)
				par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i")

				color  = densCols(num_SNP[tmp_gene], tmp_microevo_corrected_new[tmp_gene], colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
				plot(	as.numeric(num_SNP[tmp_gene]), as.numeric(tmp_microevo_corrected_new[tmp_gene]),
						xlim	= c(0, (max(num_SNP[tmp_gene]+1))), 
						xlab		= "Number of substitutions in open regions",
						ylab		= "Microevolutionary diversity (RM correction, new)",
						pch			= 16,
						col			= color,
						cex			= 0.6,
						main		= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
						las			= 1
					)

				mtext(paste(length(tmp_gene), " genes", sep=""), side = 1, line = 6, adj = 1)
				cortest_list	<- cor.test(num_SNP[tmp_gene], tmp_microevo_corrected_new[tmp_gene], method = "spearman")
				mtext(paste("microevo-Number of SNP rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)

				dev.off()
			}

		}
	}

	get_SNPInfo_PlotFigures(10000, 10000)
	get_SNPInfo_PlotFigures(5000, 5000)
	get_SNPInfo_PlotFigures(3000, 3000)
	get_SNPInfo_PlotFigures(100, 50)
}



##############################################################################################################
#### Correlation with TATA element numbers in cis-regulatory region
##############################################################################################################

if(1){
	outputdir_tmp2	<- paste(outputdir_tmp, "TATA_enrichment/", sep="")
	dir.create(outputdir_tmp2)

	stages		<- c("early", "mid", "late1", "late2")
	stage_num	<- c("15", "24", "28", "40")
	
	# Use Bed format data containing information on TATA elements detected in the upper 100 bp to lower 50 bp range of TSS.
	get_TATA_PlotFigures	<- function(upper, lower){
		outputdir_tmp3	<- paste(outputdir_tmp2, "upper", upper, "-", "lower", lower, "/", sep="")
		dir.create(outputdir_tmp3)

		for(i in 1:4){
			stage	<- stages[i]
			
			# Obtain the number of TATA elements for each gene.
			if(1){
				TATA_bed	<- read.table("/Volumes/HDD_BackUp/Fluctuation/Data/TATA_mapped_medaka-GenomeVer.2/TATAelement_closestTSS_NAomit.bed", header=FALSE, sep="\t")
				num_TATA				<- table(as.character(TATA_bed[,8]))

				output_path_tmp = paste(outputdir_tmp3, "Num-TATA_", stage, ".text", sep="")
				write.table(num_TATA, file = output_path_tmp, row.names=FALSE, col.names=FALSE, quote=FALSE)
			}


			# Visualization (vs developmental stability)
			if(1){
				tmp_flu_corrected			<- get(paste("flu_RMcorrection_", stage, sep=""))
				tmp_gene					<- intersect(names(num_TATA), names(tmp_flu_corrected))

				tmp_matrix									<- matrix(nrow=length(tmp_flu_corrected), ncol=2)
				rownames(tmp_matrix)						<- names(tmp_flu_corrected)
				tmp_matrix[,1]								<- rep(0, nrow(tmp_matrix))
				tmp_matrix[names(num_TATA[tmp_gene]), 1]	<- num_TATA[tmp_gene]
				tmp_matrix[,2]								<-  tmp_flu_corrected


				outputpath	<- paste(outputdir_tmp3, "fluctuation-NumTATA_", stage, ".pdf", sep="")
				pdf(file = outputpath,	width = 7,		height = 7)
				par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1), mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i")

				color  = densCols(tmp_matrix[,1], tmp_matrix[,2], colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
				plot(	tmp_matrix[,1], tmp_matrix[,2],
						xlim	= c(-0.5, max(tmp_matrix[,1])+0.5), 
						xlab		= "Number of TATA element",
						ylab		= "Fluctuation (RM correction)",
						pch			= 16,
						col			= color,
						cex			= 0.6,
						main		= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
						las			= 1
					)

				mtext(paste(nrow(tmp_matrix), " genes", sep=""), side = 1, line = 6, adj = 1)
				cortest_list	<- cor.test(tmp_matrix[,1], tmp_matrix[,2], method = "spearman")
				mtext(paste("fluctuation-Number of TATA rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)

				dev.off()
			}

			# Visualization (vs intraspecies diversity)
			if(1){
				tmp_microevo_corrected_new	<- get(paste("microevo_RMcorrection_new_", stage, sep=""))
				tmp_gene					<- intersect(names(num_TATA), names(tmp_microevo_corrected_new))

				tmp_matrix									<- matrix(nrow=length(tmp_microevo_corrected_new), ncol=2)
				rownames(tmp_matrix)						<- names(tmp_microevo_corrected_new)
				tmp_matrix[,1]								<- rep(0, nrow(tmp_matrix))
				tmp_matrix[names(num_TATA[tmp_gene]), 1]	<- num_TATA[tmp_gene]
				tmp_matrix[,2]								<- tmp_microevo_corrected_new


				outputpath	<- paste(outputdir_tmp3, "microevo-NumTATA_", stage, ".pdf", sep="")
				pdf(file = outputpath,	width = 7,		height = 7)
				par(xpd = TRUE, mar = c(8, 6, 4.1, 2.1),　mgp = c(5, 2, 0), pty="s", xaxs="i", yaxs="i")

				color  = densCols(tmp_matrix[,1], tmp_matrix[,2], colramp = colorRampPalette(c("#0288D1","#03A9F4","#4FC3F7","#81D4FA")), nbin = 3000)
				plot(	tmp_matrix[,1], tmp_matrix[,2],
						xlim	= c(-0.5, max(tmp_matrix[,1])+0.5), 
						xlab		= "Number of TATA element",
						ylab		= "Microevolutionary diversity (RM correction, new)",
						pch			= 16,
						col			= color,
						cex			= 0.6,
						main		= paste("After gene selection, Log, meanTPM=>",TPMCutOff, ", ", stage, sep=""),
						las			= 1
					)

				mtext(paste(nrow(tmp_matrix), " genes", sep=""), side = 1, line = 6, adj = 1)
				cortest_list	<- cor.test(tmp_matrix[,1], tmp_matrix[,2], method = "spearman")
				mtext(paste("microevo-Number of TATA rho = ", signif(cortest_list$estimate, 2), "(P = ", signif(cortest_list$p.value, 2),")", sep=""),  side = 1, line = 7, adj = 1)

				dev.off()
			}

		}
	}

	get_TATA_PlotFigures(100, 50)
}



##############################################################################################################
#### Phylogenetic tree of the biological/technical replicates
##############################################################################################################


if(1){
	outputdir_tmp2	<- paste(outputdir_tmp, "Distance_test/", sep="")
	dir.create(outputdir_tmp2)


		TPM_tmp				<- read.table("/Volumes/HDD_BackUp/Fluctuation/Data/TPM_20M_bio-tech_NOmt.txt", header=TRUE, sep="\t", row.names=1)
		TPM_tmp				<- TPM_tmp[,-1]
		TPM_adult			<- read.table("~/Documents/work/Research/Fluctuation/Data/ExpressionTables/temporal_pleiotropy/Medaka_TPM_TimeCourse.txt", header=TRUE, sep="\t", row.names=1)
		TPM_bind			<- cbind(TPM_tmp[intersect(rownames(TPM_tmp),rownames(TPM_adult)),], TPM_adult[intersect(rownames(TPM_tmp),rownames(TPM_adult)),52])
		TPM_bind			<- log(TPM_bind+1, base=10)
		colnames(TPM_bind)	<- c(	"st.23.5_bio1","st.23.5_bio2","st.23.5_bio3","st.23.5_bio4",
									"st.23.5_tech1","st.23.5_tech2","st.23.5_tech3","st.23.5_tech4",
									"st.28_tech1","st.28_tech2","st.28_tech3","st.28_tech4",
									"adult_F1"
									)
		rm(TPM_tmp)
		rm(TPM_adult)
		gc()
		gc()

	# Variance of differential gene expression
	if(1){

		tmp_dist			<- matrix(nrow=ncol(TPM_bind), ncol=ncol(TPM_bind))
		rownames(tmp_dist)	<- colnames(TPM_bind)
		colnames(tmp_dist)	<- colnames(TPM_bind)

		get_expdifference	<- function(TPM_pair){
			TPM_pair	<- unlist(TPM_pair)
			if(TPM_pair[1] == 0 && TPM_pair[2] == 0){
				return(NA)
			}else{
				if(TPM_pair[1] < TPMCutOff || TPM_pair[2] < TPMCutOff){
					return(NA)
				}else{
					return((TPM_pair[1] - TPM_pair[2]))
				}
			}
		}

		for(i in 1:ncol(TPM_bind)){
			for(k in i:ncol(TPM_bind)){
				tmp_dist[k,i]	<- var(apply(TPM_bind[,c(i,k)], 1, get_expdifference), na.rm=TRUE)		
			}
		}

		tmp_dist	<- as.dist(tmp_dist)

		outputpath	<- paste(outputdir_tmp2, "Var_AllGenes.pdf", sep="")
		pdf(file = outputpath,	width = 10,		height = 10)
	   	par(xpd = TRUE, mar = c(6, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", mfrow=c(2,2))

		tmp_cluster 		<- hclust(tmp_dist, "average")
		plot(tmp_cluster, main	="Var, All, Group-average method", las=1, xlab="", ylab="", cex.main=0.6)
		
		tmp_cluster 		<- hclust(tmp_dist, "single")
		plot(tmp_cluster, main	="Var, All, Single-linkage method", las=1, xlab="", ylab="", cex.main=0.6)

		tmp_cluster 		<- hclust(tmp_dist, "complete")
		plot(tmp_cluster, main	="Var, All, Complete-linkage method", las=1, xlab="", ylab="", cex.main=0.6)

		tmp_cluster 		<- hclust(tmp_dist, "ward.D")
		plot(tmp_cluster, main	="Var, All, Ward's method", las=1, xlab="", ylab="", cex.main=0.6)

		dev.off()
	}

	#  1- Spearman's rho
	if(1){
		tmp_dist 	<- as.dist(1 - cor(TPM_bind, method="spearman"))
	
		outputpath	<- paste(outputdir_tmp2, "1-Spearman_AllGenes", ".pdf", sep="")
		pdf(file = outputpath,	width = 10,		height = 10)
	   	par(xpd = TRUE, mar = c(6, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", mfrow=c(2,2))

		tmp_cluster 		<- hclust(tmp_dist, "average")
		plot(tmp_cluster, main	="1-Spearman, All, Group-average method", las=1, xlab="", ylab="", cex.main=0.6)
		
		tmp_cluster 		<- hclust(tmp_dist, "single")
		plot(tmp_cluster, main	="1-Spearman, All, Single-linkage method", las=1, xlab="", ylab="", cex.main=0.6)

		tmp_cluster 		<- hclust(tmp_dist, "complete")
		plot(tmp_cluster, main	="1-Spearman, All, Complete-linkage method", las=1, xlab="", ylab="", cex.main=0.6)

		tmp_cluster 		<- hclust(tmp_dist, "ward.D")
		plot(tmp_cluster, main	="1-Spearman, All, Ward's method", las=1, xlab="", ylab="", cex.main=0.6)

		dev.off()
	}



	# Bootstrap values
	if(0){
		function_var_distance	<- function(TPM_bind_Table){
			tmp_dist			<- matrix(nrow=ncol(TPM_bind), ncol=ncol(TPM_bind))
			rownames(tmp_dist)	<- colnames(TPM_bind)
			colnames(tmp_dist)	<- colnames(TPM_bind)

			get_expdifference	<- function(TPM_pair){
				TPM_pair	<- unlist(TPM_pair)
				if(TPM_pair[1] == 0 && TPM_pair[2] == 0){
					return(NA)
				}else{
					if(TPM_pair[1] < TPMCutOff || TPM_pair[2] < TPMCutOff){
						return(NA)
					}else{
						return((TPM_pair[1] - TPM_pair[2]))
					}
				}
			}	
	
			for(i in 1:ncol(TPM_bind_Table)){
				for(k in i:ncol(TPM_bind_Table)){
					tmp_dist[k,i]	<- var(apply(TPM_bind_Table[,c(i,k)], 1, get_expdifference), na.rm=TRUE)		
				}
			}

			return(as.dist(tmp_dist))
		}


		outputpath	<- paste(outputdir_tmp2, "Var_AllGenes_bootstrap,nboot=500", ".pdf", sep="")
		pdf(file = outputpath,	width = 18,		height = 18)
	   	par(xpd = TRUE, mar = c(6, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s", mfrow=c(2,2))

		result <- pvclust(TPM_bind, method.hclust="average", method.dist=function_var_distance, nboot=500)
		plot(result, main = "Var, All, Group-average method", las=1, xlab="", ylab="", cex.main=0.6)

		result <- pvclust(TPM_bind, method.hclust="single", method.dist=function_var_distance, nboot=500)
		plot(result, main = "Var, All, Single-linkage method", las=1, xlab="", ylab="", cex.main=0.6)

		result <- pvclust(TPM_bind, method.hclust="complete", method.dist=function_var_distance, nboot=500)
		plot(result, main = "Var, All, Complete-linkage method", las=1, xlab="", ylab="", cex.main=0.6)

		result <- pvclust(TPM_bind, method.hclust="ward.D", method.dist=function_var_distance, nboot=500)
		plot(result, main = "Var, All, Ward's method", las=1, xlab="", ylab="", cex.main=0.6)

		dev.off()
	}

}


##############################################################################################################
#### Depth dependency
##############################################################################################################

if(1){
	outputdir_tmp2	<- paste(outputdir_tmp, "DepthDependency/", sep="")
	dir.create(outputdir_tmp2)

	depth			<- c("3000000", "5000000", "10000000", "15000000", "20000000", "25000000", "30000000")

	# 
	if(1){
		biodistance		<- matrix(nrow=6, ncol=7)
		techdistance	<- matrix(nrow=6, ncol=7)
		for(i in 1:7){
			tmp_depth			<- depth[i]

			tmp_tablepath		<- paste("/Volumes/HDD_BackUp/Fluctuation/Data/Medaka_TPM_NOmt_", tmp_depth,"_replicates.txt", sep="")
			TPM_tmp				<- read.table(tmp_tablepath, header=TRUE, sep="\t", row.names=1)
			TPM_tmp				<- TPM_tmp[,-1]
			TPM_tmp				<- log(TPM_tmp+1, base=10)

			biodistance[,i]		<- Make_dist_array_tech(TPM_tmp[,1:4])
			techdistance[,i]	<- Make_dist_array_tech(TPM_tmp[,5:8])
		}

		colnames(biodistance)	<- paste(as.numeric(depth)/1000000, "", sep="")
		colnames(techdistance)	<- paste(as.numeric(depth)/1000000, "", sep="")	
	}

	if(1){
		outputpath	<- paste(outputdir_tmp2, "Depthdependency_VAR.pdf", sep="")
		pdf(file = outputpath,	width = 7,		height = 7)
	   	par(xpd = TRUE, mar = c(6, 4, 4, 2),  mgp = c(3.5, 1, 0), pty="s")
	
		box_position	<- boxplot(	biodistance[,1], biodistance[,2], biodistance[,3], biodistance[,4], biodistance[,5], biodistance[,6], biodistance[,7],
									col		= "deepskyblue",
									border	= "grey35",
					    			names	= colnames(biodistance),
					    			ylim	= c(0, 0.04),
									xlab	= "Read depth (M)",  
									ylab	= "VAR of TPM difference",
									cex		= 0.3,
									cex.main= 0.8,
									pch 	= 16,
			       					main	= "Log, All genes",
									las		= 1
					       			)
										
		par(new=T)

		for(i in 1:7){
			assign(paste("dist_",i, sep=""), biodistance[,i])
		}
			
		dist				<- c(dist_1,dist_2,dist_3,dist_4,dist_5,dist_6,dist_7)
		stage				<- c(rep("3M",6),rep("5M",6),rep("10M",6),rep("15M",6),rep("20M",6),rep("25M",6),rep("30M",6))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("3M","5M","10M","15M","20M","25M","30M"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= T,
					xaxt		= "n",
					ylim		= c(0, 0.04),
					ylab		= "",
					yaxt		= "n",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(0/255,133/255,201/255, alpha=0.5)), 
					)

		# Statistical test Kruskal Wallis)
		KW_result <- kruskal.test(x=list(dist_1,dist_2,dist_3,dist_4,dist_5,dist_6,dist_7))
		mtext(paste("p-values (Kruskal-Wallis, bio replicates)    ", KW_result$p.value, sep=""), side = 1, line = 4, adj = 1, cex = 0.8)

		par(new=T)

		box_position	<- boxplot(	techdistance[,1], techdistance[,2], techdistance[,3], techdistance[,4], techdistance[,5], techdistance[,6], techdistance[,7],
									col		= "gray50",
									border	= "grey35",
					       			names	= colnames(techdistance),
					       			ylim	= c(0, 0.04),
									cex		= 0.3,
									pch 	= 16,
					       			xaxt	= "n",
					       			yaxt	= "n",
					       			main	= ""
					       		)
										
		par(new=T)

		for(i in 1:7){
			assign(paste("dist_",i, sep=""), techdistance[,i])
		}
			
		dist				<- c(dist_1,dist_2,dist_3,dist_4,dist_5,dist_6,dist_7)
		stage				<- c(rep("3M",6),rep("5M",6),rep("10M",6),rep("15M",6),rep("20M",6),rep("25M",6),rep("30M",6))
		DistStage			<- data.frame(dist,stage)
		DistStage$stage		<- factor( DistStage$stage, levels(stage)<-c("3M","5M","10M","15M","20M","25M","30M"))

		stripchart(	DistStage$dist ~ DistStage$stage,
					method 		= "jitter", 
					vertical 	= T,
					xaxt		= "n",
					ylim		= c(0, 0.04),
					ylab		= "",
					yaxt		= "n",
					cex			= 1.5,
					pch 		= 16, 
					col 		= c(bg=rgb(131/255,131/255,131/255, alpha=0.5)), 
					)
		# Statistical test Kruskal Wallis)
		KW_result <- kruskal.test(x=list(dist_1,dist_2,dist_3,dist_4,dist_5,dist_6,dist_7))
		mtext(paste("p-values (Kruskal-Wallis, tech replicates)    ", KW_result$p.value, sep=""), side = 1, line = 5, adj = 1, cex = 0.8)

		dev.off()
	}

}
