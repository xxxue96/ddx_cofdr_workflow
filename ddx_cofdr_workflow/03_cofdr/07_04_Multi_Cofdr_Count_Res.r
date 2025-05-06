library(openxlsx)
library(tidyverse)
library(janitor)
library(data.table)
source("Excel_Format_Params.r")

reformat_res_multi <- function(pheno, sheets)
{
	print(paste(pheno, "start!"))
	dat = sheets[[which(names(sheets)==pheno)]]
	
	# T1.sig_snps
	if(is.null(dat$T1.sig_snps)){
		sig_snps = 0
	} else{
		sig_snps = dat$T1.sig_snps[-1, ] %>% row_to_names(row_number = 1) %>% select(SNP)
		sig_snps = nrow(sig_snps)
	}
	
	# T2.risk_loci
	if(is.null(dat$T2.risk_loci)){
		risk_loci = 0
	} else{
		risk_loci = dat$T2.risk_loci[-1, ] %>% row_to_names(row_number = 1) %>% select(rsID)
		risk_loci = nrow(risk_loci)
	}
	
	# T7.commonLoci_cofdr_hypr
	if(is.null(dat$T7.commonLoci_cofdr_hypr)){
		commonLoci_cofdr_hypr = 0
	} else{
		commonLoci_cofdr_hypr = dat$T7.commonLoci_cofdr_hypr[-1, ] %>% row_to_names(row_number = 1) %>% filter(total_detect==2)%>% select(SNP)
		commonLoci_cofdr_hypr = nrow(commonLoci_cofdr_hypr)
	}
	
	# T8.gene_based_cofdr
	if(is.null(dat$T8.gene_based_cofdr)){
		gene_based_cofdr = mgene = egene = sgene = pgene = 0
	} else{
		gene_based_cofdr = dat$T8.gene_based_cofdr[-1, ] %>% row_to_names(row_number = 1) %>% filter(total_detect>=3) %>% select(hgnc_symbol)
		gene_based_cofdr = nrow(unique(gene_based_cofdr))
		
		mgene = dat$T8.gene_based_cofdr[-1, ] %>% row_to_names(row_number = 1) %>% filter(magma_detect==1) %>% select(hgnc_symbol)
		mgene = nrow(unique(mgene))
		
		egene = dat$T8.gene_based_cofdr[-1, ] %>% row_to_names(row_number = 1) %>% filter(twas_detect==1) %>% select(hgnc_symbol)
		egene = nrow(unique(egene))
		
		sgene = dat$T8.gene_based_cofdr[-1, ] %>% row_to_names(row_number = 1) %>% filter(sptwas_detect==1) %>% select(hgnc_symbol)
		sgene = nrow(unique(sgene))
		
		pgene = dat$T8.gene_based_cofdr[-1, ] %>% row_to_names(row_number = 1) %>% filter(pwas_detect==1) %>% select(hgnc_symbol)
		pgene = nrow(unique(pgene))
	}
	
	return(cbind(pheno, sig_snps, risk_loci, mgene, egene, sgene, pgene, gene_based_cofdr, commonLoci_cofdr_hypr))
}


Multi_Cofdr_Count_Res <- function(excel_dir, save_path)
{
	paths = list.files(excel_dir, full.names=TRUE)	
	phenos = stringi::stri_replace_all_regex(basename(paths), pattern=".xlsx", replacement="", vectorize=FALSE)
	
	sheets = pbapply::pblapply(paths, read_all_sheets)	
	names(sheets) = phenos
	res = do.call(rbind.data.frame, lapply(phenos, reformat_res_multi, sheets=sheets))
	
	fwrite(res, file=save_path, sep="\t", quote=FALSE)
}


# Multi_Cofdr_Count_Res(excel_dir=paste0("03_cofdr/11_Multi_Cofdr_Summary/Tables/", c("A2", "B2", "C2")), save_path="03_cofdr/11_Multi_Cofdr_Summary/Tables/COVID_Cofdr_Count.txt")
