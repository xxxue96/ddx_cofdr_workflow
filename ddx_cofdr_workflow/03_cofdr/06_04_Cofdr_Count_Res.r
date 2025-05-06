library(openxlsx)
library(tidyverse)
library(janitor)
library(data.table)
source("Excel_Format_Params.r")

reformat_res <- function(pheno, sheets)
{
	print(paste(pheno, "start!"))
	dat = sheets[[which(names(sheets)==pheno)]]
	
	# T1.sig_snps
	if(is.null(dat$T1.sig_snps)){
		sig_snps = 0
	} else{
		sig_snps = dat$T1.sig_snps[-1, ] %>% row_to_names(row_number = 1) %>% dplyr::select(SNP)
		sig_snps = nrow(sig_snps)
	}
	
	# T2.risk_loci
	if(is.null(dat$T2.risk_loci)){
		risk_loci = 0
	} else{
		risk_loci = dat$T2.risk_loci[-1, ] %>% row_to_names(row_number = 1) %>% dplyr::select(rsID)
		risk_loci = nrow(risk_loci)
	}
	
	# T7.commonLoci_cofdr_gwaspw_hypr
	if(is.null(dat$T7.commonLoci_cofdr_gwaspw_hypr)){
		commonLoci_cofdr_gwaspw_hypr = 0
	} else{
		commonLoci_cofdr_gwaspw_hypr = dat$T7.commonLoci_cofdr_gwaspw_hypr[-1, ] %>% row_to_names(row_number = 1) %>% filter(total_detect==3)%>% dplyr::select(SNP)
		commonLoci_cofdr_gwaspw_hypr = nrow(commonLoci_cofdr_gwaspw_hypr)
	}
	
	# T8.gene_based_cofdr
	if(is.null(dat$T8.gene_based_cofdr)){
		magma_gene = twas_gene = sptwas_gene = pwas_gene = gene_based_cofdr = 0
	} else{
		magma_gene = dat$T8.gene_based_cofdr[-1, ] %>% row_to_names(row_number = 1) %>% filter(magma_detect==1) %>% dplyr::select(hgnc_symbol)
		magma_gene = nrow(unique(magma_gene))
		
		twas_gene = dat$T8.gene_based_cofdr[-1, ] %>% row_to_names(row_number = 1) %>% filter(twas_detect==1) %>% dplyr::select(hgnc_symbol)
		twas_gene = nrow(unique(twas_gene))
		
		sptwas_gene = dat$T8.gene_based_cofdr[-1, ] %>% row_to_names(row_number = 1) %>% filter(sptwas_detect==1) %>% dplyr::select(hgnc_symbol)
		sptwas_gene = nrow(unique(sptwas_gene))
		
		pwas_gene = dat$T8.gene_based_cofdr[-1, ] %>% row_to_names(row_number = 1) %>% filter(pwas_detect==1) %>% dplyr::select(hgnc_symbol)
		pwas_gene = nrow(unique(pwas_gene))
		
		gene_based_cofdr = dat$T8.gene_based_cofdr[-1, ] %>% row_to_names(row_number = 1) %>% mutate(total_detect=as.numeric(total_detect)) %>% filter(total_detect>=3) %>% dplyr::select(hgnc_symbol)
		gene_based_cofdr = nrow(unique(gene_based_cofdr))
	}
	
	# T9.sig_gene_GSEA
	if(is.null(dat$T9.sig_gene_GSEA)){
		sig_gene_GSEA = 0
	} else{
		sig_gene_GSEA = dat$T9.sig_gene_GSEA[-1, ] %>% row_to_names(row_number = 1) %>% dplyr::select(GeneSet)
		sig_gene_GSEA = nrow(unique(sig_gene_GSEA))
	}
	
	return(cbind(pheno, sig_snps, risk_loci, magma_gene, twas_gene, sptwas_gene, pwas_gene, gene_based_cofdr, sig_gene_GSEA, commonLoci_cofdr_gwaspw_hypr))
}

Cofdr_Count_Res <- function(excel_dir, save_path)
{
	paths = list.files(excel_dir, pattern=".xlsx", full.names=TRUE)	
	phenos = stringi::stri_replace_all_regex(basename(paths), pattern=c("_0.8TT.xlsx", "_no_sigsnp.xlsx", ".xlsx"), replacement=c("", "", ""), vectorize=FALSE)
	
	sheets = pbapply::pblapply(paths, read_all_sheets, cl=40)	
	names(sheets) = phenos
	res = do.call(rbind.data.frame, lapply(phenos, reformat_res, sheets=sheets))
	
	fwrite(res, file=save_path, sep="\t", quote=FALSE)
}


# Cofdr_Count_Res(excel_dir=paste0("03_cofdr/10_Excel_Summary/Tables/", c("A2", "B2", "C2"), "/respiratory"), save_path="03_cofdr/10_Excel_Summary/Tables/COVID_Cofdr_Count.txt")
