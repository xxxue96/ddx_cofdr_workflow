# annot gwascatalog results downloaded from FUMA with TT detected by co-fdr
library(data.table)
library(doParallel)
library(foreach)
library(tidyverse)
library(dplyr)
setDTthreads(20)
registerDoParallel(5)

Annot_gwas_catalog <- function(annot_path, catalog_path, save_path)
{
	dat = fread(annot_path)
	if("TT" %in% colnames(dat)){dat = dat[, c("SNP", "TT")]} else{dat = dat[, c("SNP", "P")]}
	colnames(dat) = c("snp", "TT")
	catalog = fread(catalog_path)

	# annotate IndSigSNP
	merged = merge(catalog, dat, by.x = "IndSigSNP", by.y = "snp")
	colnames(merged)[ncol(merged)] = "TT_IndSigSNP"

	# annotate single snps
	res1 = merge(merged, dat, by = "snp")

	# annotate snps containing "; "
	res2 = merged[grep("; ", merged$snp), ]
	res2 <- foreach(i = 1:nrow(res2), .combine = rbind) %dopar% {		
	
		tmp = res2[i, ]
		tmp = as.data.frame(separate_rows(tmp, snp, sep="; "))
		tmp = list(tmp,dat) %>% reduce(left_join, by="snp")
		tmp = as.data.frame(tmp %>% group_by(across(c(-snp, -TT))) %>% summarise(across(c(snp, TT), ~ paste(.x, collapse="; "))))
		tmp
	}

	# annotate snps containing " x "
	res3 = merged[grep(" x ", merged$snp), ]
	res3 <- foreach(i = 1:nrow(res3), .combine = rbind) %dopar% {
		
		tmp = res3[i, ]
		tmp = as.data.frame(separate_rows(tmp, snp, sep=" x "))
		tmp = list(tmp,dat) %>% reduce(left_join, by="snp")
		tmp = as.data.frame(tmp %>% group_by(across(c(-snp, -TT))) %>% summarise(across(c(snp, TT), ~ paste(.x, collapse=" x "))))
		tmp
	}

	# annotate snps not in gwas summary statistics but in 1KG ref genome
	res = rbind(res1, res2, res3)
	res4 = merged[!(merged$snp %in% res$snp), ]
	res4$TT = "NA"

	# combine all annotated results
	res = rbind(res, res4)
	res = res[!is.na(res$IndSigSNP), ]
	colnames(res)[ncol(res)] = "TT_snp"

	# reformat
	# column_order = match(c("GenomicLocus", "IndSigSNP", "TT_IndSigSNP", "chr", "bp", "snp", "TT_snp", "Trait"), colnames(res))
	res = res[, c(3,2,38,4,5,1,39,13,6:12,14:37)]
	res$TT_snp = as.numeric(res$TT_snp)
	res = res[order(res$TT_snp, decreasing=TRUE), ]

	# output
	fwrite(res, file = save_path, sep = "\t", quote=FALSE, na=NA)
	return(catalog_path)
}
# exmple uasge:
# annot_path = "02_ddx/01_DDx_Results/A2_Pneumonia.meta.txt.gz"
# catalog_path = "02_ddx/04_FUMA_Results/gwascatalog.txt"
# Annot_gwas_catalog(annot_path, catalog_path, save_path=catalog_path)
