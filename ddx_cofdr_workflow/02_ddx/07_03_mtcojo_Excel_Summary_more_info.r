library(openxlsx)
library(tidyverse)
library(janitor)
library(data.table)
source("Excel_Format_Params.r")
source("02_ddx/mtCOJO_Workflow_params.r")

reformat_snp <- function(pheno, snp)
{
	if(nrow(snp)>0){
		snp = as.data.frame(cbind(pheno, snp))
		snp = snp[order(as.numeric(snp$P), decreasing=FALSE),]
	} else{
		snp[1,] = NA
		snp = cbind(pheno=NA, snp)
		snp = snp[0,]		
	}	
	## rename column 
	pos = grep("BETA_", colnames(snp))
	header = c("BETA_covid19","SE_covid19", "P_covid19")
	colnames(snp)[pos[1]:(pos[1]+length(header)-1)] = header	
	colnames(snp)[1] = "Comparisons"
	
	return(snp)
}

reformat_gene <- function(pheno, gene)
{
	if(nrow(gene)>0){
		gene = as.data.frame(pheno) %>% cbind(gene)
		gene = gene[order(as.numeric(gene$P), decreasing=FALSE),]
	} else{
		gene[1,] = NA
		gene = as.data.frame(pheno) %>% cbind(gene)
		gene = gene[0,]		
	}
	## rename column 
	pos = grep("P_", colnames(gene))
	header = c("P_covid19", "FDR_covid19", "Z_covid19")
	if(identical(grep("Z_", colnames(gene)), integer(0))){
		header = c("P_covid19", "FDR_covid19")
	}
	colnames(gene)[pos[1]:(pos[1]+length(header)-1)] = header
	colnames(gene)[1] = "Comparisons"
	
	return(gene)
}

reformat_gene <- function(pheno, gene)
{
	# magam_detect
	magma_ddx = gene %>% filter(magma_detect==1) %>% select(Gene,grep("_magma", colnames(gene), value=TRUE))
	if(nrow(magma_ddx)>0){
		magma_ddx = as.data.frame(cbind(pheno, unique(magma_ddx)))		
		magma_ddx = magma_ddx[order(as.numeric(magma_ddx$P_mtCOJO_magma), decreasing=FALSE),]
	} else{
		magma_ddx[1,] = NA
		magma_ddx = cbind(pheno=NA, magma_ddx)
		magma_ddx = magma_ddx[0,]		
	}
	colnames(magma_ddx) = gsub("_magma", "", colnames(magma_ddx))
	pos = grep("Z_", colnames(magma_ddx)); colnames(magma_ddx)[pos[2]] = "Z_covid19"
	pos = grep("P_", colnames(magma_ddx)); colnames(magma_ddx)[pos[2]] = "P_covid19"
	pos = grep("FDR_", colnames(magma_ddx)); colnames(magma_ddx)[pos[2]] = "FDR_covid19"
	colnames(magma_ddx)[1] = "Comparisons"
		
	# twas_detect
	twas_ddx = gene %>% filter(twas_detect==1) %>% select(Gene,grep("_twas", colnames(gene), value=TRUE))
	if(nrow(twas_ddx)>0){
		twas_ddx = as.data.frame(cbind(pheno, unique(twas_ddx)))							
		twas_ddx = twas_ddx[order(as.numeric(twas_ddx$P_mtCOJO_twas), decreasing=FALSE),]
	} else{
		twas_ddx[1,] = NA
		twas_ddx = cbind(pheno=NA, twas_ddx)
		twas_ddx = twas_ddx[0,]	
	}
	colnames(twas_ddx) = gsub("_twas", "", colnames(twas_ddx))
	pos = grep("Z_", colnames(twas_ddx)); colnames(twas_ddx)[pos[2]] = "Z_covid19"
	pos = grep("P_", colnames(twas_ddx)); colnames(twas_ddx)[pos[2]] = "P_covid19"
	pos = grep("FDR_", colnames(twas_ddx)); colnames(twas_ddx)[pos[2]] = "FDR_covid19"
	colnames(twas_ddx)[1] = "Comparisons"
		
	# sptwas_detect
	sptwas_ddx = gene %>% filter(sptwas_detect==1) %>% select(Gene,grep("_sptwas", colnames(gene), value=TRUE))
	if(nrow(sptwas_ddx)>0){
		sptwas_ddx = as.data.frame(cbind(pheno, unique(sptwas_ddx)))							
		sptwas_ddx = sptwas_ddx[order(as.numeric(sptwas_ddx$P_mtCOJO_sptwas), decreasing=FALSE),]
	} else{
		sptwas_ddx[1,] = NA
		sptwas_ddx = cbind(pheno=NA, sptwas_ddx)
		sptwas_ddx = sptwas_ddx[0,]	
	}
	colnames(sptwas_ddx) = gsub("_sptwas", "", colnames(sptwas_ddx))
	pos = grep("Z_", colnames(sptwas_ddx)); colnames(sptwas_ddx)[pos[2]] = "Z_covid19"
	pos = grep("P_", colnames(sptwas_ddx)); colnames(sptwas_ddx)[pos[2]] = "P_covid19"
	pos = grep("FDR_", colnames(sptwas_ddx)); colnames(sptwas_ddx)[pos[2]] = "FDR_covid19"
	colnames(sptwas_ddx)[1] = "Comparisons"
	
	#pwas_detect
	pwas_ddx = gene %>% filter(pwas_detect==1) %>% select(Gene,grep("_pwas", colnames(gene), value=TRUE))
	if(nrow(pwas_ddx)>0){
		pwas_ddx = as.data.frame(cbind(pheno, unique(pwas_ddx)))		
		pwas_ddx = pwas_ddx[order(as.numeric(pwas_ddx$P_mtCOJO_pwas), decreasing=FALSE),]
	} else{
		pwas_ddx[1,] = NA
		pwas_ddx = cbind(pheno=NA, pwas_ddx)
		pwas_ddx = pwas_ddx[0,]		
	}
	colnames(pwas_ddx) = gsub("_pwas", "", colnames(pwas_ddx))
	pos = grep("Z_", colnames(pwas_ddx)); colnames(pwas_ddx)[pos[2]] = "Z_covid19"
	pos = grep("P_", colnames(pwas_ddx)); colnames(pwas_ddx)[pos[2]] = "P_covid19"
	pos = grep("FDR_", colnames(pwas_ddx)); colnames(pwas_ddx)[pos[2]] = "FDR_covid19"
	colnames(pwas_ddx)[1] = "Comparisons"
	
	return(list(magma_ddx, twas_ddx, sptwas_ddx, pwas_ddx))
}

reformat_set <- function(pheno, set)
{
	if(nrow(set)>0){
		set = as.data.frame(cbind(pheno, set))
		set = set[order(as.numeric(set$P_mtCOJO), decreasing=FALSE),]
	} else{
		set[1,] = NA
		set = cbind(pheno=NA, set)
		set = set[0,]		
	}	
	## rename column 
	pos = grep("P_", colnames(set)); colnames(set)[pos[2]] = "P_covid19"
	pos = grep("FDR_", colnames(set)); colnames(set)[pos[2]] = "FDR_covid19"	
	colnames(set)[1] = "Comparisons"
			
	return(set)
}


format_sig_snp_more_info <- function(pheno, sheets)
{
	print(paste(pheno, "start!"))
	dat = sheets[[which(names(sheets)==pheno)]]
	
	# T1.sig_snps
	sig_snps = dat$T1.sig_snps[-1, ] %>% row_to_names(row_number = 1)
	sig_snps = reformat_snp(pheno, sig_snps)
	
	# T7.gene_mtcojo
	gene_mtcojo = dat$T7.gene_mtcojo[-1, ] %>% row_to_names(row_number = 1) %>% filter(total_detect>=3) %>% unique
	gene_mtcojo = reformat_gene(pheno, gene_mtcojo)
	
	# T8.magma_gene_set
	magma_gene_set = dat$T8.magma_gene_set[-1, ] %>% row_to_names(row_number = 1) %>% filter(FDR_mtCOJO<0.05)
	magma_gene_set = reformat_set(pheno, magma_gene_set)
	
	# T9.magma_tissue
	magma_tissue = dat$T9.magma_tissue[-1, ] %>% row_to_names(row_number = 1) %>% filter(FDR_mtCOJO<0.05)
	magma_tissue = reformat_set(pheno, magma_tissue)
	
	# T10.magma_cell_type
	magma_cell_type = dat$T10.magma_cell_type[-1, ] %>% row_to_names(row_number = 1) %>% filter(FDR_mtCOJO<0.05)
	magma_cell_type = reformat_set(pheno, magma_cell_type)
	
	tmp = c(list(sig_snps), gene_mtcojo, list(magma_gene_set), list(magma_tissue), list(magma_cell_type))
	return(tmp)
}

mtcojo_Excel_Summary_more_info <- function(excel_dir, save_path)
{	
	paths = list.files(excel_dir, pattern=".xlsx", full.names=TRUE)	
	phenos = stringi::stri_replace_all_regex(basename(paths), pattern=".xlsx", replacement="", vectorize=FALSE)
	
	sheets = lapply(paths, read_all_sheets)	
	names(sheets) = phenos
	res = lapply(phenos, format_sig_snp_more_info, sheets=sheets)
	
	wb <- createWorkbook()
	addWorksheet(wb, "T1.sig_snps")
	addWorksheet(wb, "T2.validate_gene_magma")
	addWorksheet(wb, "T3.validate_gene_twas")
	addWorksheet(wb, "T4.validate_gene_sptwas")
	addWorksheet(wb, "T5.validate_gene_pwas")
	addWorksheet(wb, "T6.sig_gene_set")
	addWorksheet(wb, "T7.sig_tissue")
	addWorksheet(wb, "T8.sig_cell_type")


	h1 = paste("Table 1.Significant SNPs detected by mtCOJO at p-value threshold", pthres, "across all multi-trait comparisons")
	subh1 = 'SNP=rsid of the significant SNPs detected by mtCOJO; A1=non-effect allele; A2=effect allele; FRQ=effect allele frequency; BETA=effect size of "SNP" estimated by mtCOJO; SE=standard error of "SE"; P=pvalue of "SNP" estimated by mtCOJO; Z=zscore of "SNP" estimated by mtCOJO; N=total sample size involved in mtCOJO analysis; Other columns refer to GWAS for each trait involved in mtCOJO analysis.'
	dat1 = as.data.frame(rbindlist(lapply(res, '[[', 1), fill=TRUE))
	dat1[, c(3:4,7:ncol(dat1))] <- sapply(dat1[, c(3:4,7:ncol(dat1))], as.numeric)
	format_header(wb, "T1.sig_snps", c(1, ncol(dat1)), h1, header_style)
	format_header(wb, "T1.sig_snps", c(2, ncol(dat1)), subh1, subheader_style)
	format_main(wb, "T1.sig_snps", dat1, title_style, main_style)
	
	h2 = "Table 2. MAGMA-mtCOJO analysis results for all validated genes across all multi-trait comparisons"
	subh2 = 'Gene=HGNC symbol of the differential gene; P=the gene p-value estimated by MAGMA based on mtCOJO-derived GWAS; FDR=adjusted P; Other columns refer to the (adjusted) gene p-value estimated by MAGMA based on GWAS of each trait involved in mtCOJO analysis. "NA" represents that the "Gene" is not identified by MAGMA based on GWAS of that trait.'
	dat2 = as.data.frame(rbindlist(lapply(res, '[[', 2), fill=TRUE))
	dat2[, -c(1:2)] <- sapply(dat2[, -c(1:2)], as.numeric)
	format_header(wb, "T2.validate_gene_magma", c(1, ncol(dat2)), h2, header_style)
	format_header(wb, "T2.validate_gene_magma", c(2, ncol(dat2)), subh2, subheader_style, heights=13*2)
	format_main(wb, "T2.validate_gene_magma", dat2, title_style, main_style)
	
	h3 = "Table 3. TWAS-mtCOJO analysis results for all validated genes across all multi-trait comparisons"
	subh3 = 'Gene=HGNC symbol of the differential gene; tissue=the eQTL model of every single tissue used to obtain gene-based z-scores in TWAS analysis; P=the gene p-value estimated by TWAS based on mtCOJO-derived GWAS; FDR=adjusted P; Z=the gene z-score estimated by TWAS based on mtCOJO-derived GWAS; Other columns refer to the (adjusted) gene p-value or z-score estimated by TWAS based on GWAS of each trait involved in mtCOJO analysis. "NA" represents that the "Gene" is not identified by TWAS based on GWAS of that trait.'
	dat3 = as.data.frame(rbindlist(lapply(res, '[[', 3), fill=TRUE))
	dat3[, -c(1:3)] <- sapply(dat3[, -c(1:3)], as.numeric)
	format_header(wb, "T3.validate_gene_twas", c(1, ncol(dat3)), h3, header_style)
	format_header(wb, "T3.validate_gene_twas", c(2, ncol(dat3)), subh3, subheader_style, heights=13*2)
	format_main(wb, "T3.validate_gene_twas", dat3, title_style, main_style)
	
	h4 = "Table 4. spTWAS-mtCOJO analysis results for all validated genes across all multi-trait comparisons"
	subh4 = 'Gene=HGNC symbol of the differential gene; tissue=the sQTL model of every single tissue used to obtain gene-based z-scores in spTWAS analysis; P=the gene p-value estimated by spTWAS based on mtCOJO-derived GWAS; FDR=adjusted P; Z=the gene z-score estimated by spTWAS based on mtCOJO-derived GWAS; Other columns refer to the (adjusted) gene p-value or z-score estimated by spTWAS based on GWAS of each trait involved in mtCOJO analysis. "NA" represents that the "Gene" is not identified by spTWAS based on GWAS of that trait.'
	dat4 = as.data.frame(rbindlist(lapply(res, '[[', 4), fill=TRUE))
	dat4[, -c(1:3)] <- sapply(dat4[, -c(1:3)], as.numeric)
	format_header(wb, "T4.validate_gene_sptwas", c(1, ncol(dat4)), h4, header_style)
	format_header(wb, "T4.validate_gene_sptwas", c(2, ncol(dat4)), subh4, subheader_style, heights=13*2)
	format_main(wb, "T4.validate_gene_sptwas", dat4, title_style, main_style)
	
	h5 = "Table 5. PWAS-mtCOJO analysis results for all validated genes across all multi-trait comparisons"
	subh5 = 'Gene=HGNC symbol of the differential gene; P=the gene p-value estimated by PWAS based on mtCOJO-derived GWAS; FDR=adjusted P; Z=the gene z-score estimated by PWAS based on mtCOJO-derived GWAS; Other columns refer to the (adjusted) gene p-value or z-score estimated by PWAS based on GWAS of each trait involved in mtCOJO analysis. "NA" represents that the "Gene" is not identified by PWAS based on GWAS of that trait.'
	dat5 = as.data.frame(rbindlist(lapply(res, '[[', 5), fill=TRUE))
	dat5[, -c(1:2)] <- sapply(dat5[, -c(1:2)], as.numeric)
	format_header(wb, "T5.validate_gene_pwas", c(1, ncol(dat5)), h5, header_style)
	format_header(wb, "T5.validate_gene_pwas", c(2, ncol(dat5)), subh5, subheader_style, heights=13*2)
	format_main(wb, "T5.validate_gene_pwas", dat5, title_style, main_style)
	
	h6 = "Table 6. Gene set analysis based on mtCOJO-derived GWAS using MAGMA"
	subh6 = 'Gene_set=full name of the differential gene set; P_mtCOJO=the gene set p-value estimated by MAGMA based on mtCOJO-derived GWAS; FDR_mtCOJO=adjusted P_mtCOJO; P_trait1=the gene set p-value estimated by MAGMA based on GWAS of trait1; FDR_trait1=adjusted P_trait1; P_trait2=the gene set p-value estimated by MAGMA based on GWAS of trait2; FDR_trait2=adjusted P_trait2'
	dat6 = as.data.frame(rbindlist(lapply(res, '[[', 6), fill=TRUE))
	dat6[, -c(1:2)] <- sapply(dat6[, -c(1:2)], as.numeric)
	format_header(wb, "T6.sig_gene_set", c(1, ncol(dat6)), h6, header_style)
	format_header(wb, "T6.sig_gene_set", c(2, ncol(dat6)), subh6, subheader_style, heights=13*4)
	format_main(wb, "T6.sig_gene_set", dat6, title_style, main_style)
	
	h7 = "Table 7. Tissue enrichment analysis based on mtCOJO-derived GWAS using MAGMA"
	subh7 = 'Tissue=full name of the differential tissue; P_mtCOJO=the tissue p-value estimated by MAGMA based on mtCOJO-derived GWAS; FDR_mtCOJO=adjusted P_mtCOJO; P_trait1=the tissue p-value estimated by MAGMA based on GWAS of trait1; FDR_trait1=adjusted P_trait1; P_trait2=the tissue p-value estimated by MAGMA based on GWAS of trait2; FDR_trait2=adjusted P_trait2'
	dat7 = as.data.frame(rbindlist(lapply(res, '[[', 7), fill=TRUE))
	dat7[, -c(1:2)] <- sapply(dat7[, -c(1:2)], as.numeric)
	format_header(wb, "T7.sig_tissue", c(1, ncol(dat7)), h7, header_style)
	format_header(wb, "T7.sig_tissue", c(2, ncol(dat7)), subh7, subheader_style, heights=13*4)
	format_main(wb, "T7.sig_tissue", dat7, title_style, main_style)
	
	h8 = "Table 8. Significant cell types based on mtCOJO-derived GWAS using MAGMA"
	subh8 = 'Dataset=scRNA-seq dataset with analyzed cell types; Cell_type=full name of the differential cell type; P_mtCOJO=the cell type p-value estimated by MAGMA based on mtCOJO-derived GWAS; FDR_mtCOJO=adjusted P_mtCOJO; P_trait1=the cell type p-value estimated by MAGMA based on GWAS of trait1; FDR_trait1=adjusted P_trait1; P_trait2=the cell type p-value estimated by MAGMA based on GWAS of trait2; FDR_trait2=adjusted P_trait2'
	dat8 = as.data.frame(rbindlist(lapply(res, '[[', 8), fill=TRUE))
	dat8[, -c(1:3)] <- sapply(dat8[, -c(1:3)], as.numeric)
	format_header(wb, "T8.sig_cell_type", c(1, ncol(dat8)), h8, header_style)
	format_header(wb, "T8.sig_cell_type", c(2, ncol(dat8)), subh8, subheader_style, heights=13*4)
	format_main(wb, "T8.sig_cell_type", dat8, title_style, main_style)

	saveWorkbook(wb, file = save_path, overwrite = TRUE)
	
	return(save_path)
}



# covidhgi = c("A2", "B2", "C2")
# t2dir = "respiratory"

# mtCOJO_Excel_Summary(excel_dir=paste0("02_ddx/13_mtCOJO_Excel_Summary/Tables/", covidhgi), save_path="02_ddx/13_mtCOJO_Excel_Summary/Tables/Summary_A2_B2_C2_respiratory_more_info.xlsx")




