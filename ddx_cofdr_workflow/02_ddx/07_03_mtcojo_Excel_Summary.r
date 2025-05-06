library(openxlsx)
library(tidyverse)
library(janitor)
source("Excel_Format_Params.r")
source("02_ddx/mtCOJO_Workflow_params.r")

format_sig_snp <- function(pheno, sheets)
{
	print(paste(pheno, "start!"))
	dat = sheets[[which(names(sheets)==pheno)]]
	
	# T1.sig_snps
	sig_snps = dat$T1.sig_snps[-1, ] %>% row_to_names(row_number = 1) %>% select(SNP)
	sig_snps = as.data.frame(cbind(sig_snps, pheno))
	colnames(sig_snps) = c("SNP", "comparisons")
	
	# T2.risk_loci
	risk_loci = dat$T2.risk_loci[-1, ] %>% row_to_names(row_number = 1) %>% select(rsID)
	risk_loci = as.data.frame(cbind(risk_loci, pheno))
	colnames(risk_loci) = c("SNP", "comparisons")
	
	# T7.gene_mtcojo
	gene_mtcojo = dat$T7.gene_mtcojo[-1, ] %>% row_to_names(row_number = 1)
	if(nrow(gene_mtcojo)>0){
		gene_mtcojo = gene_mtcojo %>% cbind(pheno) %>% filter(total_detect>=3) %>% select(Gene, pheno) %>% unique
	} else{
		gene_mtcojo = as.data.frame(matrix(ncol=2, nrow=0))		
	}	
	colnames(gene_mtcojo) = c("Gene", "comparisons")
	
	# T8.magma_gene_set
	magma_gene_set = dat$T8.magma_gene_set[-1, ] %>% row_to_names(row_number = 1) 
	if(nrow(magma_gene_set)>0){
		magma_gene_set = magma_gene_set %>% cbind(pheno) %>% filter(FDR_mtCOJO<0.05) %>% select(Gene_set, pheno) %>% unique
	} else{
		magma_gene_set = as.data.frame(matrix(ncol=2, nrow=0))		
	}
	colnames(magma_gene_set) = c("Gene_set", "comparisons")
	
	# T9.magma_tissu
	magma_tissue = dat$T9.magma_tissu[-1, ] %>% row_to_names(row_number = 1)
	if(nrow(magma_tissue)>0){
		magma_tissue = magma_tissue %>% cbind(pheno) %>% filter(FDR_mtCOJO<0.05) %>% select(Tissue, pheno) %>% unique				
	} else{
		magma_tissue = as.data.frame(matrix(ncol=2, nrow=0))		
	}
	colnames(magma_tissue) = c("Tissue", "comparisons")
	
	# T10.magma_cell_type
	magma_cell_type = dat$T10.magma_cell_type[-1, ] %>% row_to_names(row_number = 1)
	if(nrow(magma_cell_type)>0){
		magma_cell_type = magma_cell_type %>% cbind(pheno) %>% filter(FDR_mtCOJO<0.05) %>% select(Cell_type, pheno) %>% unique
	} else{
		magma_cell_type = as.data.frame(matrix(ncol=2, nrow=0))		
	}
	colnames(magma_cell_type) = c("Cell_type", "comparisons")
	
	return(list(sig_snps, risk_loci, gene_mtcojo, magma_gene_set, magma_tissue, magma_cell_type))
}

mtcojo_Excel_Summary <- function(excel_dir, save_path)
{
	paths = list.files(excel_dir, pattern=".xlsx", full.names=TRUE)	
	phenos = stringi::stri_replace_all_regex(basename(paths), pattern=".xlsx", replacement="", vectorize=FALSE)
	
	sheets = lapply(paths, read_all_sheets)	
	names(sheets) = phenos
	res = lapply(phenos, format_sig_snp, sheets=sheets)
	
	# idx: index for "sig_snps, risk_loci, gene_mtcojo, magma_gene_set, magma_tissue, magma_cell_type"
	# col1: name of the 1st column
	format_Excel <- function(idx, col1)
	{
		print(idx)
		dat = lapply(res, '[[', idx) 
		dat = dat[lapply(dat, nrow)>0] 
		
		if(!identical(dat, list())){
			dat = dat %>% reduce(full_join, by=col1)
			if(ncol(dat)==2){
				dat$N_comparisons = 1
			}else{
				dat$N_comparisons = apply(dat[, -1], 1, function(x) sum(!is.na(x)))
				dat$comparisons = apply(dat[, -c(1,ncol(dat))], 1, paste, collapse = ";")
			}
			dat = dat[, c(col1, "N_comparisons", "comparisons")]
			dat = dat[order(dat$N_comparisons, decreasing=TRUE),]
		} else{
			dat = setNames(as.data.frame(matrix(ncol=3, nrow=1)), c(col1, "N_comparisons", "comparisons"))
		}
			return(dat)
	}
	
	wb <- createWorkbook()
	addWorksheet(wb, "T1.sig_snps")
	addWorksheet(wb, "T2.risk_loci")
	addWorksheet(wb, "T3.validate_gene_mtcojo")
	addWorksheet(wb, "T4.magma_gene_set")
	addWorksheet(wb, "T5.magma_tissue")
	addWorksheet(wb, "T6.magma_cell_type")


	h1 = paste("Table 1.Significant SNPs detected by mtCOJO at p-value threshold", pthres, "across all multi-trait comparisons")
	subh1 = 'SNP=rsid of the significant SNP; N_comparisons=number of comparisons identifying the "SNP"; comparisons=name of comparisons identifying the "SNP"'
	dat1 = format_Excel(1, "SNP")
	format_header(wb, "T1.sig_snps", c(1, 13), h1, header_style)
	format_header(wb, "T1.sig_snps", c(2, 13), subh1, subheader_style)
	format_main(wb, "T1.sig_snps", dat1, title_style, main_style)
	
	h2 = "Table 2. Identification of genomic risk loci across all multi-trait comparisons"
	subh2 = 'SNP=rsid of the lead SNP in the genomic locus; N_comparisons=number of comparisons identifying the "SNP"; comparisons=name of comparisons identifying the "SNP"'
	dat2 = format_Excel(2, "SNP")
	format_header(wb, "T2.risk_loci", c(1, 13), h2, header_style)
	format_header(wb, "T2.risk_loci", c(2, 13), subh2, subheader_style)
	format_main(wb, "T2.risk_loci", dat2, title_style, main_style)
		
	h3 = "Table 3. Results summary of validated genes associated with conditional traits across all multi-trait comparisons"
	subh3 = 'Gene=name of the significant gene; N_comparisons=number of comparisons identifying the "Gene"; comparisons=name of comparisons identifying the "Gene"'
	dat3 = format_Excel(3, "Gene")
	format_header(wb, "T3.validate_gene_mtcojo", c(1, 12), h3, header_style)
	format_header(wb, "T3.validate_gene_mtcojo", c(2, 12), subh3, subheader_style)
	format_main(wb, "T3.validate_gene_mtcojo", dat3, title_style, main_style)
	
	h4 = "Table 4.  Gene set analysis based on mtCOJO-derived GWAS using MAGMA across all multi-trait comparisons"
	subh4 = 'Gene_set=name of the significant gene set; N_comparisons=number of comparisons identifying the "Gene_set"; comparisons=name of comparisons identifying the "Gene_set"'
	dat4 = format_Excel(4, "Gene_set")
	format_header(wb, "T4.magma_gene_set", c(1, 12), h4, header_style)
	format_header(wb, "T4.magma_gene_set", c(2, 12), subh4, subheader_style)
	format_main(wb, "T4.magma_gene_set", dat4, title_style, main_style)
	
	h5 = "Table 5. Tissue enrichment analysis based on mtCOJO-derived GWAS using MAGMA across all multi-trait comparisons"
	subh5 = 'Tissue=name of the significant tissue; N_comparisons=number of comparisons identifying the "Tissue"; comparisons=name of comparisons identifying the "Tissue"'
	dat5 = format_Excel(5, "Tissue")
	format_header(wb, "T5.magma_tissue", c(1, 13), h5, header_style)
	format_header(wb, "T5.magma_tissue", c(2, 13), subh5, subheader_style)
	format_main(wb, "T5.magma_tissue", dat5, title_style, main_style)
	
	h6 = "Table 6.  Significant cell types based on mtCOJO-derived GWAS using MAGMA per dataset cell type analysis across all multi-trait comparisons"
	subh6 = 'Cell_type=name of the significant cell type; N_comparisons=number of comparisons identifying the "Cell_type"; comparisons=name of comparisons identifying the "Cell_type"'
	dat6 = format_Excel(6, "Cell_type")
	format_header(wb, "T6.magma_cell_type", c(1, 14), h6, header_style)
	format_header(wb, "T6.magma_cell_type", c(2, 14), subh6, subheader_style)
	format_main(wb, "T6.magma_cell_type", dat6, title_style, main_style)
	
	saveWorkbook(wb, file = save_path, overwrite = TRUE)
	
	return(save_path)
}
# covidhgi = c("A2", "B2", "C2")
# t2phenos = c("Asthma", "ILD", "IPF", "COPD", "Pneumonia.meta", "Pneumonia.b.meta", "Pneumonia.v.meta", "Flu.meta")

# mtcojo_Excel_Summary(excel_dir=paste0("02_ddx/13_mtCOJO_Excel_Summary/Tables/", covidhgi), save_path="02_ddx/13_mtCOJO_Excel_Summary/Tables/Summary_A2_B2_C2_respiratory.xlsx")


