library(openxlsx)
library(tidyverse)
library(janitor)
library(data.table)
setDTthreads(40)
source("Excel_Format_Params.r")
source("02_ddx/mtCOJO_Workflow_params.r")

reformat_res_mtcojo <- function(pheno, sheets)
{
	print(paste(pheno, "start!"))
	dat = sheets[[which(names(sheets)==pheno)]]
	
	# T1.sig_snps
	sig_snps = dat$T1.sig_snps[-1, ] %>% row_to_names(row_number = 1) %>% select(SNP)
	sig_snps = nrow(sig_snps)
	
	# T2.risk_loci
	risk_loci = dat$T2.risk_loci[-1, ] %>% row_to_names(row_number = 1) %>% select(rsID)
	risk_loci = nrow(risk_loci)
	
	# T7.gene_mtcojo
	magma_ddx = dat$T7.gene_mtcojo[-1, ] %>% row_to_names(row_number = 1) %>% filter(magma_detect==1) %>% select(Gene)
	magma_ddx = nrow(unique(magma_ddx))
		
	twas_ddx = dat$T7.gene_mtcojo[-1, ] %>% row_to_names(row_number = 1) %>% filter(twas_detect==1) %>% select(Gene)
	twas_ddx = nrow(unique(twas_ddx))
		
	sptwas_ddx = dat$T7.gene_mtcojo[-1, ] %>% row_to_names(row_number = 1) %>% filter(sptwas_detect==1) %>% select(Gene)
	sptwas_ddx = nrow(unique(sptwas_ddx))
		
	pwas_ddx = dat$T7.gene_mtcojo[-1, ] %>% row_to_names(row_number = 1) %>% filter(pwas_detect==1) %>% select(Gene)
	pwas_ddx = nrow(unique(pwas_ddx))
	
	validate_gene = dat$T7.gene_mtcojo[-1, ] %>% row_to_names(row_number = 1) %>% filter(total_detect>=3) %>% select(Gene)
	validate_gene = nrow(unique(validate_gene))
	
	# T8.magma_gene_set
	magma_gene_set = dat$T8.magma_gene_set[-1, ] %>% row_to_names(row_number = 1) %>% filter(FDR_mtCOJO<FDRthres2) %>% select(Gene_set)
	magma_gene_set = nrow(unique(magma_gene_set))
	
	# T9.magma_tissue
	magma_tissue = dat$T9.magma_tissue[-1, ] %>% row_to_names(row_number = 1) %>% filter(FDR_mtCOJO<FDRthres2) %>% select(Tissue)
	magma_tissue = nrow(magma_tissue)
	
	# T10.magma_cell_type
	magma_cell_type = dat$T10.magma_cell_type[-1, ] %>% row_to_names(row_number = 1) %>% filter(FDR_mtCOJO<FDRthres2) %>% select(Cell_type)
	magma_cell_type = nrow(unique(magma_cell_type))
	
	unadj_snps = fread(paste0("01_data/00_Standardise_GWAS/",t1,".txt.gz")) %>% filter(P<pthres) %>% nrow
	
	tmp = dat$T1.sig_snps[-1, ] %>% row_to_names(row_number = 1) %>% mutate(P = as.numeric(P))
	tmp[, paste0("P_",t1)] = as.numeric(tmp[, paste0("P_",t1)])	
	specific_snps = tmp[tmp[, paste0("P_",t1)]<pthres & tmp$P<tmp[, paste0("P_",t1)],] %>% nrow
		
	return(cbind(pheno, unadj_snps, sig_snps, risk_loci, specific_snps, magma_ddx, twas_ddx, sptwas_ddx, pwas_ddx, validate_gene, magma_gene_set, magma_tissue, magma_cell_type))
}

mtcojo_Count_Res <- function(excel_dir, save_path)
{
	paths = list.files(excel_dir, pattern=".xlsx", full.names=TRUE)	
	phenos = stringi::stri_replace_all_regex(basename(paths), pattern=".xlsx", replacement="", vectorize=FALSE)
	
	sheets = pbapply::pblapply(paths, read_all_sheets, cl=40)	
	names(sheets) = phenos
	res = do.call(rbind.data.frame, lapply(phenos, reformat_res_mtcojo, sheets=sheets))
	
	fwrite(res, file=save_path, sep="\t", quote=FALSE)
	return(save_path)
}

# mtcojo_Count_Res(excel_dir="02_ddx/13_mtCOJO_Excel_Summary/Tables", save_path="02_ddx/13_mtCOJO_Excel_Summary/Tables/COVID_mtcojo_Count.txt")
