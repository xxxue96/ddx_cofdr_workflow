library(tidyverse)
library(dplyr)
library(data.table)


annot_magma <- function(traits, magma)
{
	ann = lapply(traits, function(trait){
							dat = fread(paste0("01_data/03_Gene_Result/",trait, "/", trait, ".0UP.0DOWN_sig.magma.txt.gz"))[, c("hgnc_symbol", "P", "FDR", "Z")]
							colnames(dat)[2:4] = paste0(c("P_", "FDR_", "Z_"), trait)
							return(dat)})
	magma = c(list(magma), ann) %>% reduce(inner_join, by="hgnc_symbol")
	return(magma)
}

annot_twas <- function(traits, twas)
{
	ann = lapply(traits, function(trait){
							dat = fread(paste0("01_data/03_Gene_Result/",trait, "/", trait, "_sig.eqtl.txt.gz"))[, c("hgnc_symbol", "tissue", "pvalue", "FDR", "zscore")]
							colnames(dat)[3:5] = paste0(c("P_", "FDR_", "Z_"), trait)
							return(dat)})
	## only keep the gene with largest TT
	twas = c(list(twas), ann) %>% reduce(inner_join, by=c("hgnc_symbol", "tissue")) %>% group_by(hgnc_symbol, tissue) %>% slice(which.max(TT)) %>% as.data.table %>% unique
	return(twas)
}

annot_sptwas <- function(traits, sptwas)
{
	ann = lapply(traits, function(trait){
							dat = fread(paste0("01_data/03_Gene_Result/",trait, "/", trait, "_sig.sqtl.txt.gz"))[, c("hgnc_symbol", "tissue", "pvalue", "FDR", "zscore", "sqtl_hg38")]
							colnames(dat)[3:5] = paste0(c("P_", "FDR_", "Z_"), trait)
							return(dat)})
	sptwas = c(list(sptwas), ann) %>% reduce(inner_join, by=c("hgnc_symbol", "tissue", "sqtl_hg38")) %>% dplyr::select(-sqtl_hg38) %>% group_by(hgnc_symbol, tissue) %>% slice(which.max(TT)) %>% as.data.table %>% unique
	return(sptwas)
}

annot_pwas <- function(traits, pwas)
{
	ann = lapply(traits, function(trait){
							dat = fread(paste0("01_data/03_Gene_Result/",trait, "/", trait, "_sig.pqtl.txt.gz"))[, c("hgnc_symbol", "PWAS.P", "FDR", "PWAS.Z")]
							colnames(dat)[2:4] = paste0(c("P_", "FDR_", "Z_"), trait)
							return(dat)})
	pwas = c(list(pwas), ann) %>% reduce(inner_join, by="hgnc_symbol") %>% group_by(hgnc_symbol) %>% slice(which.max(TT)) %>% as.data.table %>% unique
	return(pwas)
}

# compare genes from magma-,twas-,sptwas-,and pwas-based cofdr results
# annotate identified genes with magma-,twas-,sptwas-,and pwas-results of original GWAS for analyzed traits
Gene_Compare <- function(traits, magma_path, twas_path, sptwas_path, pwas_path, save_path)
{
	magma_gene = fread(magma_path)[, c("hgnc_symbol", "TT")] %>% mutate_at(1, as.character) %>% annot_magma(traits = traits)	
	eqtl_gene = fread(twas_path)[, c("hgnc_symbol", "TT", "tissue")] %>% mutate_at(c(1,3), as.character) %>% annot_twas(traits = traits)
	sqtl_gene = fread(sptwas_path)[, c("hgnc_symbol", "TT", "tissue", "sqtl_hg38")] %>% mutate_at(c(1,3,4), as.character) %>% annot_sptwas(traits = traits)
	pqtl_gene = fread(pwas_path)[, c("hgnc_symbol", "TT")] %>% mutate_at(1, as.character) %>% annot_pwas(traits = traits)
	
	# merge
	colnames(magma_gene)[-1] = paste0(colnames(magma_gene)[-1], "_magma")
	colnames(eqtl_gene)[-1] = paste0(colnames(eqtl_gene)[-1], "_twas")
	colnames(sqtl_gene)[-1] = paste0(colnames(sqtl_gene)[-1], "_sptwas")
	colnames(pqtl_gene)[-1] = paste0(colnames(pqtl_gene)[-1], "_pwas")
	all_gene = list(magma_gene, eqtl_gene, sqtl_gene, pqtl_gene) %>% reduce(full_join, by = "hgnc_symbol")
	
	# compare
	if(nrow(all_gene)>0)
	{
		magma_detect = twas_detect = sptwas_detect = pwas_detect = total_detect = numeric(nrow(all_gene))	
		attach(all_gene)
		for(i in 1:nrow(all_gene))
		{
			if(!is.na(TT_magma[i])) {magma_detect[i] = 1} else{magma_detect[i] = 0}
			if(!is.na(TT_twas[i])) {twas_detect[i] = 1} else{twas_detect[i] = 0}
			if(!is.na(TT_sptwas[i])) {sptwas_detect[i] = 1} else{sptwas_detect[i] = 0}
			if(!is.na(TT_pwas[i])) {pwas_detect[i] = 1} else{pwas_detect[i] = 0}

			total_detect[i] = sum(magma_detect[i], twas_detect[i], sptwas_detect[i], pwas_detect[i])
		}
		detach(all_gene)	
		all_gene = cbind(all_gene, magma_detect, twas_detect, sptwas_detect, pwas_detect, total_detect)	
		
		
		idx = c(as.character(outer(paste0("Z_",traits), c("_magma","_twas","_sptwas","_pwas"), paste0)), as.character(outer(paste0("P_",traits), c("_magma","_twas","_sptwas","_pwas"), paste0)), as.character(outer(paste0("FDR_",traits), c("_magma","_twas","_sptwas","_pwas"), paste0)))
		all_gene = cbind(all_gene[, -idx, with=FALSE], all_gene[, idx, with=FALSE])
		all_gene = all_gene[order(all_gene$total_detect, decreasing=TRUE),]
		
		# save
		fwrite(all_gene, file=save_path, sep="\t", na=NA, quote=FALSE)
	} else{
		print("no magma/eqtl/sqtl/pqtl gene exists!")
	}
}

Run_Gene_Compare <- function(t1, t2, res_dir)
{
	traits = c(t1, t2)
	pheno =  paste(traits, collapse="_")
	
	magma_path = paste0(res_dir,"/",pheno,".0UP.0DOWN_sig.magma.txt.gz")
	twas_path = paste0(res_dir,"/",pheno,"_sig.eqtl.txt.gz")
	sptwas_path = paste0(res_dir,"/",pheno,"_sig.sqtl.txt.gz")
	pwas_path = paste0(res_dir,"/",pheno,"_sig.pqtl.txt.gz")
	save_path = paste0(res_dir,"/",pheno,"_geneCompare.txt.gz")

	Gene_Compare(traits, magma_path, twas_path, sptwas_path, pwas_path, save_path)
}

# covidhgi = c("A2", "B2", "C2")
# respiratory = c("Asthma", "ILD", "IPF", "COPD", "Pneumonia.meta", "Pneumonia.b.meta", "Pneumonia.v.meta", "Flu.meta")
# for(t1 in covidhgi){
	# for(t2 in respiratory){
		# print(paste(t1, t2, "start"))
		# Run_Gene_Compare(t1, t2, res_dir=paste0("03_cofdr/04_Cofdr_Gene_Result/",t1,"_",t2))
	# }
# }
