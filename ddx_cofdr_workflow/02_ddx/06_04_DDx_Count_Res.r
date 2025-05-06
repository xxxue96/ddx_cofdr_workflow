library(openxlsx)
library(tidyverse)
library(janitor)
library(data.table)
source("Excel_Format_Params.r")
source("02_ddx/DDx_Workflow_params.r")

reformat_res <- function(pheno, sheets)
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
	
	# T7.commonLoci_ddx_ccgwas
	if(is.null(dat$T7.commonLoci_ddx_ccgwas)){
		commonLoci_ddx_ccgwas = 0
	} else{
		commonLoci_ddx_ccgwas = dat$T7.commonLoci_ddx_ccgwas[-1, ] %>% row_to_names(row_number = 1) %>% select(SNP)
		commonLoci_ddx_ccgwas = nrow(commonLoci_ddx_ccgwas)
	}
	
	# T8.gene_ddx
	if(is.null(dat$T8.gene_ddx)){
		magma_ddx = twas_ddx = sptwas_ddx = pwas_ddx = 0
	} else{
		magma_ddx = dat$T8.gene_ddx[-1, ] %>% row_to_names(row_number = 1) %>% filter(magma_detect==1) %>% select(Gene)
		magma_ddx = nrow(unique(magma_ddx))
		
		twas_ddx = dat$T8.gene_ddx[-1, ] %>% row_to_names(row_number = 1) %>% filter(twas_detect==1) %>% select(Gene)
		twas_ddx = nrow(unique(twas_ddx))
		
		sptwas_ddx = dat$T8.gene_ddx[-1, ] %>% row_to_names(row_number = 1) %>% filter(sptwas_detect==1) %>% select(Gene)
		sptwas_ddx = nrow(unique(sptwas_ddx))
		
		pwas_ddx = dat$T8.gene_ddx[-1, ] %>% row_to_names(row_number = 1) %>% filter(pwas_detect==1) %>% select(Gene)
		pwas_ddx = nrow(unique(pwas_ddx))
	}
	
	# T10.commonGene_ddx_ccgwas
	if(is.null(dat$T10.commonGene_ddx_ccgwas)){
		magma_gene = twas_gene = sptwas_gene = pwas_gene = commonGene_ddx_ccgwas = 0
	} else{
		magma_gene = dat$T10.commonGene_ddx_ccgwas[-1, ] %>% row_to_names(row_number = 1) %>% filter(magma_detect==1) %>% select(Gene)
		magma_gene = nrow(unique(magma_gene))
		
		twas_gene = dat$T10.commonGene_ddx_ccgwas[-1, ] %>% row_to_names(row_number = 1) %>% filter(twas_detect==1) %>% select(Gene)
		twas_gene = nrow(unique(twas_gene))
		
		sptwas_gene = dat$T10.commonGene_ddx_ccgwas[-1, ] %>% row_to_names(row_number = 1) %>% filter(sptwas_detect==1) %>% select(Gene)
		sptwas_gene = nrow(unique(sptwas_gene))
		
		pwas_gene = dat$T10.commonGene_ddx_ccgwas[-1, ] %>% row_to_names(row_number = 1) %>% filter(pwas_detect==1) %>% select(Gene)
		pwas_gene = nrow(unique(pwas_gene))
		
		commonGene_ddx_ccgwas = dat$T10.commonGene_ddx_ccgwas[-1, ] %>% row_to_names(row_number = 1) %>% filter(total_detect>=3) %>% select(Gene)
		commonGene_ddx_ccgwas = nrow(unique(commonGene_ddx_ccgwas))
	}
		
	# T11.magma_gene_set
	magma_gene_set = dat$T11.magma_gene_set[-1, ] %>% row_to_names(row_number = 1) %>% filter(FDR_DDx<FDRthres2) %>% select(Gene_set)
	magma_gene_set = nrow(unique(magma_gene_set))
	
	# T12.magma_tissue
	magma_tissue = dat$T12.magma_tissue[-1, ] %>% row_to_names(row_number = 1) %>% filter(FDR_DDx<FDRthres2) %>% select(Tissue)
	magma_tissue = nrow(magma_tissue)
	
	# T13.magma_cell_type
	magma_cell_type = dat$T13.magma_cell_type[-1, ] %>% row_to_names(row_number = 1) %>% filter(FDR_DDx<FDRthres2) %>% select(Cell_type)
	magma_cell_type = nrow(unique(magma_cell_type))
	
	return(cbind(pheno, sig_snps, risk_loci, magma_ddx, twas_ddx, sptwas_ddx, pwas_ddx, magma_gene_set, magma_tissue, magma_cell_type, commonLoci_ddx_ccgwas, magma_gene, twas_gene, sptwas_gene, pwas_gene, commonGene_ddx_ccgwas))
}


DDx_Count_Res <- function(excel_dir, save_path)
{
	paths = list.files(excel_dir, pattern=".xlsx", full.names=TRUE)	
	phenos = stringi::stri_replace_all_regex(basename(paths), pattern=".xlsx", replacement="", vectorize=FALSE)
	
	sheets = pbapply::pblapply(paths, read_all_sheets, cl=40)	
	names(sheets) = phenos
	res = do.call(rbind.data.frame, lapply(phenos, reformat_res, sheets=sheets))
	
	fwrite(res, file=save_path, sep="\t", quote=FALSE)
	return(save_path)
}

# DDx_Count_Res(excel_dir=paste0("02_ddx/11_Excel_Summary/Tables/", c("A2", "B2", "C2"), "/respiratory"), save_path="02_ddx/11_Excel_Summary/Tables/COVID_DDx_Count.txt")