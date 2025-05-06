# summarize SNPs detected by different methods
library(tidyverse)
library(dplyr)
library(data.table)

# @snp_path: significantly colocalized snps detected by cofdr, gwsa-pw and hyprcoloc
# @annot_sumstats: merged sumstats of colocalized traits
SNP_Compare <- function(snp_path, annot_sumstats, save_path)
{
	# read 
	cofdr_snp = fread(snp_path[[1]])[, c("SNP", "TT", "P")]
	gwaspw_snp = fread(snp_path[[2]])[, c("SNP", "chunkPPA_3", "snpPPA_3")]
	hyprcoloc_snp = fread(snp_path[[3]])[, c("candidate_snp", "posterior_prob", "posterior_explained_by_snp")]
	colnames(hyprcoloc_snp)[1] = "SNP"
	
	# merge
	cofdr_snp$SNP=as.character(cofdr_snp$SNP)
	gwaspw_snp$SNP = as.character(gwaspw_snp$SNP)
	hyprcoloc_snp$SNP = as.character(hyprcoloc_snp$SNP)	
	all_snp = list(cofdr_snp, gwaspw_snp, hyprcoloc_snp) %>% reduce(full_join, by = "SNP")
	
	# compare
	cofdr_detect = gwaspw_detect = hyprcoloc_detect = total_detect = numeric(nrow(all_snp))	
	attach(all_snp)
	for(i in 1:nrow(all_snp))
	{
		if(!is.na(TT[i]) & !is.na(P[i])) {cofdr_detect[i] = 1} else{cofdr_detect[i] = 0}
		if(!is.na(chunkPPA_3[i]) & !is.na(snpPPA_3[i])) {gwaspw_detect[i] = 1} else{gwaspw_detect[i] = 0}
		if(!is.na(posterior_prob[i]) & !is.na(posterior_explained_by_snp[i])) {hyprcoloc_detect[i] = 1} else{hyprcoloc_detect[i] = 0}
		total_detect[i] = sum(cofdr_detect[i], gwaspw_detect[i], hyprcoloc_detect[i])
	}
	detach(all_snp)	
	all_snp = cbind(all_snp, cofdr_detect, gwaspw_detect, hyprcoloc_detect, total_detect)	
		
	# annotate
	annot_sumstats = annot_sumstats %>% relocate(sort(names(.)[-c(1:5)]))
	all_snp = merge(all_snp, annot_sumstats, by="SNP")
	all_snp = all_snp[order(all_snp$total_detect, decreasing=TRUE), ]
	
	# save
	fwrite(all_snp, file=save_path, sep="\t", na=NA, quote=FALSE)
	
	return(save_path)	
}