library(openxlsx)
library(tidyverse)
library(janitor)
library(data.table)
source("03_cofdr/06_03_Cofdr_Excel_Summary_more_info.r")

Multi_format_sig_snp <- function(pheno, sheets)
{
	print(paste(pheno, "start!"))
	dat = sheets[[which(names(sheets)==pheno)]]
	
	# T2.risk_loci
	rsID = dat$T2.risk_loci[-1, ] %>% row_to_names(row_number = 1) %>% select(rsID)
	risk_loci = dat$T7.commonLoci_cofdr_hypr[-1, ] %>% row_to_names(row_number = 1) %>% filter(SNP %in% rsID$rsID) %>% select(SNP:posterior_explained_by_snp_hyprcoloc)	
	risk_loci = reformat_snp(pheno, risk_loci)
	
	# T8.gene_based_cofdr
	if(is.null(dat$T8.gene_based_cofdr)){
		gene_based_cofdr = list(data.frame(), data.frame(), data.frame(), data.frame())
	} else{
		gene_based_cofdr = dat$T8.gene_based_cofdr[-1, ] %>% row_to_names(row_number = 1)
		gene_based_cofdr = reformat_gene(pheno, gene_based_cofdr)
	}
	
	tmp = c(list(risk_loci), gene_based_cofdr)
	return(tmp)
}

Multi_Cofdr_Excel_Summary_more_info <- function(excel_dir, save_path)
{		
	paths = list.files(excel_dir, full.names=TRUE)	
	phenos = stringi::stri_replace_all_regex(basename(paths), pattern=".xlsx", replacement="", vectorize=FALSE)
	
	sheets = pbapply::pblapply(paths, read_all_sheets, cl=40)	
	names(sheets) = phenos
	res = lapply(phenos, Multi_format_sig_snp, sheets=sheets)

	wb <- createWorkbook()
	addWorksheet(wb, "T1.risk_loci")
	addWorksheet(wb, "T2.Sig.Gene.MAGMA")
	addWorksheet(wb, "T3.Sig.Gene.TWAS")
	addWorksheet(wb, "T4.Sig.Gene.spTWAS")
	addWorksheet(wb, "T5.Sig.Gene.PWAS")
	
	h1 = "Table 1. Genomic risk loci from cofdr across all multiple trait comparisons"
	subh1 = 'A1=non-effect allele; A2=effect allele; FRQ=effect allele frequency; BETA=effect size; SE=standard error of “BETA”; P=p-value; Z=z-scores;TT_cofdr=posterior possibility from co-fdr that the corresponding SNP in the "SNP" column is colocalized; chunkPPA_3_gwaspw=posteior possibility from gwas-pw that a defined LD-independent region contains a colocalized SNP; snpPPA_3_gwaspw=posterior possibility from gwas-pw that the corresponding SNP located in the defined genomic region is the colocalized SNP; posterior_prob_hyprcoloc=posteior possibility from hyprcoloc that a defined LD-independent region contains a colocalized SNP; posterior_explained_by_snp_hyprcoloc=posterior possibility from hyprcoloc that the corresponding SNP located in the defined genomic region is the colocalized SNP;prop=the overall proportion of SNPs that are shared among traits; LRT=the p-value of the likelihood-ratio test indicating whether there is evidence of genetic overlap (i.e., shared SNPs) more than expected by chance'
	dat1 = as.data.frame(rbindlist(lapply(res, '[[', 1), fill=TRUE))
	dat1[, c(3:4,7:ncol(dat1))] <- sapply(dat1[, c(3:4,7:ncol(dat1))], as.numeric)
	format_header(wb, "T1.risk_loci", c(1, ncol(dat1)), h1, header_style)
	format_header(wb, "T1.risk_loci", c(2, ncol(dat1)), subh1, subheader_style, heights=2*13)
	format_main(wb, "T1.risk_loci", dat1, title_style, main_style)
	
	h2 = "Table 2. Genes identified by MAGMA-cofdr analysis across multi-trait colocalization analyses (Sig.Gene.MAGMA)"
	subh2 = 'Comparisons=trait comparisons identifying the "Gene"; Gene=HGNC symbol of the colocalized gene; TT=posterior possibility that corresponding Gene is shared among traits; Z=z-scores of the Gene estimated by MAGMA; prop=the overall proportion of genes that are shared among traits; LRT=the p-value of the likelihood-ratio test indicating whether there is evidence of genetic overlap (i.e., shared genes) more than expected by chance.'	
	dat2 = rbindlist(lapply(res, '[[', 2), fill=TRUE)
	dat2 = cbind(dat2[,1:2], as.data.table(apply(dat2[,-c(1:2)], 2, as.numeric)))	
	format_header(wb, "T2.Sig.Gene.MAGMA", c(1, 14), h2, header_style)
	format_header(wb, "T2.Sig.Gene.MAGMA", c(2, 14), subh2, subheader_style, heights=2*13)
	format_main(wb, "T2.Sig.Gene.MAGMA", dat2, title_style, main_style)
	
	h3 = "Table 3. Genes identified by TWAS-cofdr analysis across multi-trait colocalization analyses (Sig.Gene.TWAS)"
	subh3 = 'Comparisons=trait comparisons identifying the "Gene"; Gene=HGNC symbol of the colocalized gene; tissue=the prediction model trained on the expression data of every single tissue in GTEx v8 ("smultixcan" refers to gene-based z-scores from MultiXcan analyses across multi-tissue models); TT=posterior possibility that corresponding Gene is shared among traits; Z=z-scores of the Gene estimated by TWAS using PrediXcan; prop=the overall proportion of genes that are shared among traits; LRT=the p-value of the likelihood-ratio test indicating whether there is evidence of genetic overlap (i.e., shared genes) more than expected by chance.'	
	dat3 = rbindlist(lapply(res, '[[', 3), fill=TRUE)
	dat3 = cbind(dat3[,c(1,2,4)], as.data.table(apply(dat3[,-c(1,2,4)], 2, as.numeric)))
	format_header(wb, "T3.Sig.Gene.TWAS", c(1, 14), h3, header_style)
	format_header(wb, "T3.Sig.Gene.TWAS", c(2, 14), subh3, subheader_style, heights=3*13)
	format_main(wb, "T3.Sig.Gene.TWAS", dat3, title_style, main_style)
	
	h4 = "Table 4. Genes identified by spTWAS-cofdr analysis across multi-trait colocalization analyses (Sig.Gene.spTWAS)"
	subh4 = 'Comparisons=trait comparisons identifying the "Gene"; Gene=HGNC symbol of the colocalized gene; tissue=the prediction model trained on the splicing data of every single tissue in GTEx v8("smultixcan" refers to gene-based z-scores from MultiXcan analyses across multi-tissue models); TT=posterior possibility that corresponding Gene is shared among traits; Z=z-scores of the Gene estimated by spTWAS using PrediXcan; prop=the overall proportion of genes that are shared among traits; LRT=the p-value of the likelihood-ratio test indicating whether there is evidence of genetic overlap (i.e., shared genes) more than expected by chance.'	
	dat4 = rbindlist(lapply(res, '[[', 4), fill=TRUE)
	dat4 = cbind(dat4[,c(1,2,4)], as.data.table(apply(dat4[,-c(1,2,4)], 2, as.numeric)))	
	format_header(wb, "T4.Sig.Gene.spTWAS", c(1, 14), h4, header_style)
	format_header(wb, "T4.Sig.Gene.spTWAS", c(2, 14), subh4, subheader_style, heights=3*13)
	format_main(wb, "T4.Sig.Gene.spTWAS", dat4, title_style, main_style)
	
	h5 = "Table 5. Genes identified by PWAS-cofdr analysis across multi-trait colocalization analyses (Sig.Gene.PWAS)"
	subh5 = 'Comparisons=trait comparisons identifying the "Gene"; Gene=HGNC symbol of the colocalized gene; TT=posterior possibility that corresponding Gene is shared among traits; Z_trait1=z-scores of the Gene estimated by PWAS; prop=the overall proportion of genes that are shared among traits; LRT=the p-value of the likelihood-ratio test indicating whether there is evidence of genetic overlap (i.e., shared genes) more than expected by chance.'	
	dat5 = rbindlist(lapply(res, '[[', 5), fill=TRUE)
	dat5 = cbind(dat5[,1:2], as.data.table(apply(dat5[,-c(1:2)], 2, as.numeric)))	
	format_header(wb, "T5.Sig.Gene.PWAS", c(1, 14), h5, header_style)
	format_header(wb, "T5.Sig.Gene.PWAS", c(2, 14), subh5, subheader_style, heights=2*13)
	format_main(wb, "T5.Sig.Gene.PWAS", dat5, title_style, main_style)
	
	saveWorkbook(wb, file = save_path, overwrite = TRUE)
	
	return(save_path)
}

# Multi_Cofdr_Excel_Summary_more_info(excel_dir=paste0("03_cofdr/11_Multi_Cofdr_Summary/Tables/", c("A2", "B2", "C2")), save_path=paste0("03_cofdr/11_Multi_Cofdr_Summary/Tables/A2_B2_C2_respiratory_more_info.xlsx"))