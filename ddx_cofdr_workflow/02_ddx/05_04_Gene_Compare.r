library(tidyverse)
library(dplyr)
library(data.table)

read_sumstats <- function(path)
{
	if(is.data.frame(path)){ return(path) } else{ return(fread(path, fill = TRUE)) }
}

merge_magma <- function(ddx_res, ccgwas_res)
{
	dat1 = fread(ddx_res) %>% rename_at(vars(-hgnc_symbol), ~ paste0(., '_DDx')) %>% mutate_at("hgnc_symbol", as.character)	
	dat2 = fread(ccgwas_res) %>% rename_at(vars(-hgnc_symbol), ~ paste0(., '_OLS')) %>% mutate_at("hgnc_symbol", as.character)	
	
	return(merge(dat1, dat2, by="hgnc_symbol"))
}

merge_twas <- function(ddx_res, ccgwas_res)
{
	dat1 = fread(ddx_res) %>% rename(P=pvalue, Z=zscore) %>% rename_at(vars(-hgnc_symbol, -tissue), ~ paste0(., '_DDx'))%>% mutate_at(c("hgnc_symbol","tissue"), as.character) 	
	dat2 = fread(ccgwas_res) %>% rename(P=pvalue, Z=zscore) %>% rename_at(vars(-hgnc_symbol, -tissue), ~ paste0(., '_OLS'))%>% mutate_at(c("hgnc_symbol","tissue"), as.character) 

	return(merge(dat1, dat2, by=c("hgnc_symbol", "tissue")))
}

merge_sptwas <- function(ddx_res, ccgwas_res)
{
	dat1 = fread(ddx_res) %>% rename(P=pvalue, Z=zscore) %>% dplyr::select(-ensemble_id) %>% rename_at(vars(-hgnc_symbol, -tissue, -sqtl_hg38, -sqtl_hg19), ~ paste0(., '_DDx'))%>% mutate_at(c("hgnc_symbol","tissue","sqtl_hg38","sqtl_hg19"), as.character) 	
	dat2 = fread(ccgwas_res) %>% rename(P=pvalue, Z=zscore) %>% dplyr::select(-ensemble_id) %>% rename_at(vars(-hgnc_symbol, -tissue, -sqtl_hg38, -sqtl_hg19), ~ paste0(., '_OLS'))%>% mutate_at(c("hgnc_symbol","tissue","sqtl_hg38","sqtl_hg19"), as.character) 

	return(merge(dat1, dat2, by=c("hgnc_symbol", "tissue", "sqtl_hg38", "sqtl_hg19")))
}

merge_pwas <- function(ddx_res, ccgwas_res)
{
	dat1 = fread(ddx_res) %>% rename(P=PWAS.P, Z=PWAS.Z) %>% dplyr::select(-FILE) %>% rename_at(vars(-hgnc_symbol), ~ paste0(., '_DDx')) %>% mutate_at("hgnc_symbol", as.character)	
	dat2 = fread(ccgwas_res) %>% rename(P=PWAS.P, Z=PWAS.Z) %>% dplyr::select(-FILE) %>% rename_at(vars(-hgnc_symbol), ~ paste0(., '_OLS')) %>% mutate_at("hgnc_symbol", as.character)

	return(merge(dat1, dat2, by="hgnc_symbol"))
}

annot_magma <- function(traits, magma)
{
	ann = lapply(traits, function(trait){
							dat = fread(paste0("01_data/03_Gene_Result/",trait, "/", trait, ".0UP.0DOWN_sig.magma.txt.gz"))[, c("hgnc_symbol", "P", "FDR", "Z")]
							colnames(dat)[2:4] = paste0(c("P_", "FDR_", "Z_"), trait)
							return(dat)})
	magma = c(list(magma), ann) %>% reduce(left_join, by="hgnc_symbol")
	return(magma)
}

annot_twas <- function(traits, twas)
{
	ann = lapply(traits, function(trait){
							dat = fread(paste0("01_data/03_Gene_Result/",trait, "/", trait, "_sig.eqtl.txt.gz"))[, c("hgnc_symbol", "tissue", "pvalue", "FDR", "zscore")]
							colnames(dat)[3:5] = paste0(c("P_", "FDR_", "Z_"), trait)
							return(dat)})
	## only keep the gene with largest TT
	twas = c(list(twas), ann) %>% reduce(left_join, by=c("hgnc_symbol", "tissue")) %>% group_by(hgnc_symbol, tissue) %>% slice(which.min(FDR_DDx)) %>% as.data.table %>% unique
	return(twas)
}

annot_sptwas <- function(traits, sptwas)
{
	ann = lapply(traits, function(trait){
							dat = fread(paste0("01_data/03_Gene_Result/",trait, "/", trait, "_sig.sqtl.txt.gz"))[, c("hgnc_symbol", "tissue", "pvalue", "FDR", "zscore", "sqtl_hg38", "sqtl_hg19")]
							colnames(dat)[3:5] = paste0(c("P_", "FDR_", "Z_"), trait)
							return(dat)})
	sptwas = c(list(sptwas), ann) %>% reduce(left_join, by=c("hgnc_symbol", "tissue", "sqtl_hg38", "sqtl_hg19")) %>% dplyr::select(-sqtl_hg38, -sqtl_hg19) %>% group_by(hgnc_symbol, tissue) %>% slice(which.min(FDR_DDx)) %>% as.data.table %>% unique
	return(sptwas)
}

annot_pwas <- function(traits, pwas)
{
	ann = lapply(traits, function(trait){
							dat = fread(paste0("01_data/03_Gene_Result/",trait, "/", trait, "_sig.pqtl.txt.gz"))[, c("hgnc_symbol", "PWAS.P", "FDR", "PWAS.Z")]
							colnames(dat)[2:4] = paste0(c("P_", "FDR_", "Z_"), trait)
							return(dat)})
	pwas = c(list(pwas), ann) %>% reduce(left_join, by="hgnc_symbol") %>% group_by(hgnc_symbol) %>% slice(which.min(FDR_DDx)) %>% as.data.table %>% unique
	return(pwas)
}

# compare genes from magma-,twas-,sptwas-,and pwas-based cofdr results
# annotate identified genes with magma-,twas-,sptwas-,and pwas-results of original GWAS for analyzed traits
Gene_Compare <- function(traits, magma_path, twas_path, sptwas_path, pwas_path, save_path)
{
	magma_gene = read_sumstats(magma_path) %>% annot_magma(traits = traits)	
	eqtl_gene = read_sumstats(twas_path) %>% annot_twas(traits = traits)
	sqtl_gene = read_sumstats(sptwas_path) %>% annot_sptwas(traits = traits)
	pqtl_gene = read_sumstats(pwas_path) %>% annot_pwas(traits = traits)
	
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
			if(!is.na(P_DDx_magma[i])) {magma_detect[i] = 1} else{magma_detect[i] = 0}
			if(!is.na(P_DDx_twas[i])) {twas_detect[i] = 1} else{twas_detect[i] = 0}
			if(!is.na(P_DDx_sptwas[i])) {sptwas_detect[i] = 1} else{sptwas_detect[i] = 0}
			if(!is.na(P_DDx_pwas[i])) {pwas_detect[i] = 1} else{pwas_detect[i] = 0}

			total_detect[i] = sum(magma_detect[i], twas_detect[i], sptwas_detect[i], pwas_detect[i])
		}
		detach(all_gene)	
		all_gene = cbind(all_gene, magma_detect, twas_detect, sptwas_detect, pwas_detect, total_detect)			
	} else{
		print("no magma/eqtl/sqtl/pqtl gene exists!")
		all_gene = cbind(all_gene, data.frame(magma_detect=NA, twas_detect=NA, sptwas_detect=NA, pwas_detect=NA, total_detect=NA))
	}
	
	idx = c(as.character(outer(paste0("Z_",traits), c("_magma","_twas","_sptwas","_pwas"), paste0)), as.character(outer(paste0("P_",traits), c("_magma","_twas","_sptwas","_pwas"), paste0)), as.character(outer(paste0("FDR_",traits), c("_magma","_twas","_sptwas","_pwas"), paste0)))
	all_gene = cbind(all_gene[, -idx, with=FALSE], all_gene[, idx, with=FALSE])
	all_gene = all_gene[order(all_gene$total_detect, decreasing=TRUE),]
		
	# save
	fwrite(all_gene, file=save_path, sep="\t", na=NA, quote=FALSE)
	return(save_path)
}

 