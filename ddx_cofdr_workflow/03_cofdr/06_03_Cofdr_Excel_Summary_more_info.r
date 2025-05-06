library(openxlsx)
library(tidyverse)
library(janitor)
library(data.table)
source("Excel_Format_Params.r")

read_log_file <- function(path)
{
	log_file = fread(path)[, c("piTT", "LRTp")]
	trait = gsub(".log", "", basename(path))
	log_file = cbind(trait, log_file)
	
	return(log_file)
}

reformat_snp <- function(pheno, snp)
{
	if(nrow(snp)>0){
		snp = as.data.frame(cbind(pheno, snp))
		snp = snp[order(as.numeric(snp$TT_cofdr), decreasing=TRUE),]
		snp = cbind(snp, fread(paste0("03_cofdr/01_Cofdr_Result/",pheno,".log"))[, c("piTT", "LRTp")])
	} else{
		snp[1,] = NA
		snp = cbind(pheno=NA, snp, piTT=NA, LRTp=NA)
		snp = snp[0,]		
	}	
	## rename column 
	pos = grep("FRQ_", colnames(snp))
	header = as.vector(outer(c("FRQ_", "BETA_","SE_", "P_", "Z_"), paste0("trait", 1:length(pos)), paste0))
	colnames(snp)[pos[1]:(pos[1]+length(header)-1)] = header
	
	colnames(snp)[1] = "Comparisons"
	colnames(snp)[(length(colnames(snp))-1):length(colnames(snp))] = c("prop", "LRT")
		
	return(snp)
}

reformat_gene <- function(pheno, gene)
{
	# magam_detect
	magma_cofdr = gene %>% filter(magma_detect==1) %>% dplyr::select(hgnc_symbol,grep("_magma", colnames(gene), value=TRUE))
	if(nrow(magma_cofdr)>0){
		magma_cofdr = as.data.frame(cbind(pheno, unique(magma_cofdr)))		
		magma_cofdr = magma_cofdr[order(as.numeric(magma_cofdr$TT_magma), decreasing=TRUE),]
		magma_cofdr = cbind(magma_cofdr, fread(paste0("03_cofdr/04_Cofdr_Gene_Result/",pheno,"/01_MAGMA/",pheno,".0UP.0DOWN.log"))[, c("piTT", "LRTp")])		
	} else{
		magma_cofdr[1,] = NA
		magma_cofdr = cbind(pheno=NA, magma_cofdr, piTT=NA, LRTp=NA)
		magma_cofdr = magma_cofdr[0,]		
	}
	pos = grep("Z_", colnames(magma_cofdr)); colnames(magma_cofdr)[pos] = paste0("Z_trait", 1:length(pos))
	pos = grep("P_", colnames(magma_cofdr)); colnames(magma_cofdr)[pos] = paste0("P_trait", 1:length(pos))	
	pos = grep("FDR_", colnames(magma_cofdr)); colnames(magma_cofdr)[pos] = paste0("FDR_trait", 1:length(pos))		
	colnames(magma_cofdr)[1:3] = c("Comparisons", "Gene", "TT")
	colnames(magma_cofdr)[(length(colnames(magma_cofdr))-1):length(colnames(magma_cofdr))] = c("prop", "LRT")
		
	# twas_detect
	twas_cofdr = gene %>% filter(twas_detect==1) %>% dplyr::select(hgnc_symbol,grep("_twas", colnames(gene), value=TRUE))
	if(nrow(twas_cofdr)>0){
		twas_cofdr = as.data.frame(cbind(pheno, unique(twas_cofdr)))				

		paths = c(list.files(paste0("03_cofdr/04_Cofdr_Gene_Result/",pheno,"/02_spredixcan/eqtl"), pattern=".log", full.names=TRUE), paste0("03_cofdr/04_Cofdr_Gene_Result/",pheno,"/03_smultixcan/eqtl/",pheno,"_smultixcan.log"))
		log_file = do.call(rbind.data.frame, lapply(paths, read_log_file))
		twas_cofdr$trait = paste0(twas_cofdr$pheno, "_", twas_cofdr$tissue_twas)
		twas_cofdr = merge(twas_cofdr, log_file, by="trait")[, -1]
			
		twas_cofdr = twas_cofdr[order(as.numeric(twas_cofdr$TT_twas), decreasing=TRUE),]
	} else{
		twas_cofdr[1,] = NA
		twas_cofdr = cbind(pheno=NA, twas_cofdr, piTT=NA, LRTp=NA)
		twas_cofdr = twas_cofdr[0,]	
	}
	pos = grep("Z_", colnames(twas_cofdr)); colnames(twas_cofdr)[pos] = paste0("Z_trait", 1:length(pos))
	pos = grep("P_", colnames(twas_cofdr)); colnames(twas_cofdr)[pos] = paste0("P_trait", 1:length(pos))	
	pos = grep("FDR_", colnames(twas_cofdr)); colnames(twas_cofdr)[pos] = paste0("FDR_trait", 1:length(pos))		
	colnames(twas_cofdr)[1:4] =  c("Comparisons", "Gene", "TT", "tissue")
	colnames(twas_cofdr)[(length(colnames(twas_cofdr))-1):length(colnames(twas_cofdr))] = c("prop", "LRT")
	
	# sptwas_detect
	sptwas_cofdr = gene %>% filter(sptwas_detect==1) %>% dplyr::select(hgnc_symbol,grep("_sptwas", colnames(gene), value=TRUE))
	if(nrow(sptwas_cofdr)>0){
		sptwas_cofdr = as.data.frame(cbind(pheno, unique(sptwas_cofdr)))
			
		paths = c(list.files(paste0("03_cofdr/04_Cofdr_Gene_Result/",pheno,"/02_spredixcan/sqtl"), pattern=".log", full.names=TRUE), paste0("03_cofdr/04_Cofdr_Gene_Result/",pheno,"/03_smultixcan/sqtl/",pheno,"_smultixcan.log"))
		log_file = do.call(rbind.data.frame, lapply(paths, read_log_file))
		sptwas_cofdr$trait = paste0(sptwas_cofdr$pheno, "_", sptwas_cofdr$tissue_sptwas)
		sptwas_cofdr = merge(sptwas_cofdr, log_file, by="trait")[, -1]
			
		sptwas_cofdr = sptwas_cofdr[order(as.numeric(sptwas_cofdr$TT_sptwas), decreasing=TRUE),]			
	} else{
		sptwas_cofdr[1,] = NA
		sptwas_cofdr = cbind(pheno=NA, sptwas_cofdr, piTT=NA, LRTp=NA)
		sptwas_cofdr = sptwas_cofdr[0,]	
	}
	pos = grep("Z_", colnames(sptwas_cofdr)); colnames(sptwas_cofdr)[pos] = paste0("Z_trait", 1:length(pos))
	pos = grep("P_", colnames(sptwas_cofdr)); colnames(sptwas_cofdr)[pos] = paste0("P_trait", 1:length(pos))	
	pos = grep("FDR_", colnames(sptwas_cofdr)); colnames(sptwas_cofdr)[pos] = paste0("FDR_trait", 1:length(pos))		
	colnames(sptwas_cofdr)[1:4] =  c("Comparisons", "Gene", "TT", "tissue")
	colnames(sptwas_cofdr)[(length(colnames(sptwas_cofdr))-1):length(colnames(sptwas_cofdr))] = c("prop", "LRT")
	
	#pwas_detect
	pwas_cofdr = gene %>% filter(pwas_detect==1) %>% dplyr::select(hgnc_symbol,grep("_pwas", colnames(gene), value=TRUE))
	if(nrow(pwas_cofdr)>0){
		pwas_cofdr = as.data.frame(cbind(pheno, unique(pwas_cofdr)))
		pwas_cofdr = pwas_cofdr[order(as.numeric(pwas_cofdr$TT_pwas), decreasing=TRUE),]
		pwas_cofdr = cbind(pwas_cofdr, fread(paste0("03_cofdr/04_Cofdr_Gene_Result/",pheno,"/04_pwas/",pheno,".log"))[, c("piTT", "LRTp")])			
	} else{
		pwas_cofdr[1,] = NA
		pwas_cofdr = cbind(pheno=NA, pwas_cofdr, piTT=NA, LRTp=NA)
		pwas_cofdr = pwas_cofdr[0,]	
	}
	pos = grep("Z_", colnames(pwas_cofdr)); colnames(pwas_cofdr)[pos] = paste0("Z_trait", 1:length(pos))
	pos = grep("P_", colnames(pwas_cofdr)); colnames(pwas_cofdr)[pos] = paste0("P_trait", 1:length(pos))	
	pos = grep("FDR_", colnames(pwas_cofdr)); colnames(pwas_cofdr)[pos] = paste0("FDR_trait", 1:length(pos))		
	colnames(pwas_cofdr)[1:3] = c("Comparisons", "Gene", "TT")
	colnames(pwas_cofdr)[(length(colnames(pwas_cofdr))-1):length(colnames(pwas_cofdr))] = c("prop", "LRT")
	
	return(list(magma_cofdr, twas_cofdr, sptwas_cofdr, pwas_cofdr))
}

# @sheets: all sheets from the excel generated by 06_01_Cofdr_Format_Excel.r
format_sig_snp_more_info <- function(pheno, sheets)
{
	print(paste(pheno, "start!"))
	dat = sheets[[which(names(sheets)==pheno)]]
	
	# T7.commonLoci_cofdr_gwaspw_hypr
	if(is.null(dat$T7.commonLoci_cofdr_gwaspw_hypr)){
		commonLoci_cofdr_gwaspw_hypr = data.frame()
	} else{
		commonLoci_cofdr_gwaspw_hypr = dat$T7.commonLoci_cofdr_gwaspw_hypr[-1, ] %>% row_to_names(row_number = 1) %>% filter(total_detect==3) %>% dplyr::select(SNP:posterior_explained_by_snp_hyprcoloc)
		commonLoci_cofdr_gwaspw_hypr = reformat_snp(pheno, commonLoci_cofdr_gwaspw_hypr)
	}
	
	# T8.gene_based_cofdr
	if(is.null(dat$T8.gene_based_cofdr)){
		gene_based_cofdr = list(data.frame(), data.frame(), data.frame(), data.frame())
	} else{
		gene_based_cofdr = dat$T8.gene_based_cofdr[-1, ] %>% row_to_names(row_number = 1) %>% mutate(total_detect=as.numeric(total_detect)) %>% filter(total_detect>=3)
		gene_based_cofdr = reformat_gene(pheno, gene_based_cofdr)
	}
	
	tmp = c(list(commonLoci_cofdr_gwaspw_hypr), gene_based_cofdr)
	return(tmp)
}

Cofdr_Excel_Summary_more_info <- function(excel_dir, save_path)
{		
	paths = list.files(excel_dir, pattern=".xlsx", full.names=TRUE)	
	phenos = stringi::stri_replace_all_regex(basename(paths), pattern=".xlsx", replacement="", vectorize=FALSE)
		
	sheets = pbapply::pblapply(paths, read_all_sheets, cl=40)	
	names(sheets) = phenos
	res = lapply(phenos, format_sig_snp_more_info, sheets=sheets)

	wb <- createWorkbook()
	addWorksheet(wb, "T1.commonLoci_cofdr_gwaspw_hypr")
	addWorksheet(wb, "T2.magma_cofdr")
	addWorksheet(wb, "T3.twas_cofdr")
	addWorksheet(wb, "T4.sptwas_cofdr")
	addWorksheet(wb, "T5.pwas_cofdr")
	
	h1 = "Table 1. Validated loci across pairwise colocalization analyses"
	subh1 = 'A1=non-effect allele; A2=effect allele; FRQ_trait1=effect allele frequency for trait1; BETA_trait1=effect size for trait1; SE_trait1=standard error of “BETA” for trait1; P_trait1=p-value for trait1; Z_trait1=z-scores for trait1; FRQ_trait2=effect allele frequency for trait2; BETA_trait2=effect size for trait2; SE_trait2=standard error of “BETA” for trait2; P_trait2=p-value for trait2; Z_trait2=z-scores for trait2; TT_cofdr=posterior possibility from co-fdr that the corresponding SNP in the "SNP" column is colocalized; chunkPPA_3_gwaspw=posteior possibility from gwas-pw that a defined LD-independent region contains a colocalized SNP; snpPPA_3_gwaspw=posterior possibility from gwas-pw that the corresponding SNP located in the defined genomic region is the colocalized SNP; posterior_prob_hyprcoloc=posteior possibility from hyprcoloc that a defined LD-independent region contains a colocalized SNP; posterior_explained_by_snp_hyprcoloc=posterior possibility from hyprcoloc that the corresponding SNP located in the defined genomic region is the colocalized SNP;prop=the overall proportion of SNPs that are shared across both traits; LRT=the p-value of the likelihood-ratio test indicating whether there is evidence of genetic overlap (i.e., shared SNPs) more than expected by chance'
	dat1 = as.data.frame(rbindlist(lapply(res, '[[', 1), fill=TRUE))
	dat1[, c(3:4,7:ncol(dat1))] <- sapply(dat1[, c(3:4,7:ncol(dat1))], as.numeric)
	format_header(wb, "T1.commonLoci_cofdr_gwaspw_hypr", c(1, ncol(dat1)), h1, header_style)
	format_header(wb, "T1.commonLoci_cofdr_gwaspw_hypr", c(2, ncol(dat1)), subh1, subheader_style, heights=4*13)
	format_main(wb, "T1.commonLoci_cofdr_gwaspw_hypr", dat1, title_style, main_style)
	
	h2 = "Table 2. MAGMA-cofdr analysis results for all validated genes across pairwise colocalization analyses"
	subh2 = 'Comparisons=trait comparisons identifying the "Gene"; Gene=HGNC symbol of the colocalized gene; TT=posterior possibility that corresponding Gene is shared between trait1 and trait2; Z_trait1=z-scores of the Gene estimated by MAGMA for trait1; Z_trait2=z-scores of the Gene estimated by MAGMA for trait2; prop=the overall proportion of genes that are shared across both traits; LRT=the p-value of the likelihood-ratio test indicating whether there is evidence of genetic overlap (i.e., shared genes) more than expected by chance.'	
	dat2 = rbindlist(lapply(res, '[[', 2), fill=TRUE)
	dat2 = cbind(dat2[,1:2], as.data.table(apply(dat2[,-c(1:2)], 2, as.numeric)))
	format_header(wb, "T2.magma_cofdr", c(1, 14), h2, header_style)
	format_header(wb, "T2.magma_cofdr", c(2, 14), subh2, subheader_style, heights=2*13)
	format_main(wb, "T2.magma_cofdr", dat2, title_style, main_style)
	
	h3 = "Table 3. TWAS-cofdr analysis results for all validated genes across pairwise colocalization analyses"
	subh3 = 'Comparisons=trait comparisons identifying the "Gene"; Gene=HGNC symbol of the colocalized gene; tissue=the prediction model trained on the expression data of every single tissue in GTEx v8 ("smultixcan" refers to gene-based z-scores from MultiXcan analyses across multi-tissue models); TT=posterior possibility that corresponding Gene is shared between trait1 and trait2; Z_trait1=z-scores of the Gene estimated by TWAS using PrediXcan for trait1; Z_trait2=z-scores of the Gene estimated by TWAS using PrediXcan for trait2; prop=the overall proportion of genes that are shared across both traits; LRT=the p-value of the likelihood-ratio test indicating whether there is evidence of genetic overlap (i.e., shared genes) more than expected by chance.'	
	dat3 = rbindlist(lapply(res, '[[', 3), fill=TRUE)
	dat3 = cbind(dat3[,c(1,2,4)], as.data.table(apply(dat3[,-c(1,2,4)], 2, as.numeric)))	
	format_header(wb, "T3.twas_cofdr", c(1, 14), h3, header_style)
	format_header(wb, "T3.twas_cofdr", c(2, 14), subh3, subheader_style, heights=3*13)
	format_main(wb, "T3.twas_cofdr", dat3, title_style, main_style)
	
	h4 = "Table 4. spTWAS-cofdr analysis results for all validated genes across pairwise colocalization analyses"
	subh4 = 'Comparisons=trait comparisons identifying the "Gene"; Gene=HGNC symbol of the colocalized gene; tissue=the prediction model trained on the splicing data of every single tissue in GTEx v8("smultixcan" refers to gene-based z-scores from MultiXcan analyses across multi-tissue models); TT=posterior possibility that corresponding Gene is shared between trait1 and trait2; Z_trait1=z-scores of the Gene estimated by spTWAS using PrediXcan for trait1; Z_trait2=z-scores of the Gene estimated by spTWAS using PrediXcan for trait2; prop=the overall proportion of genes that are shared across both traits; LRT=the p-value of the likelihood-ratio test indicating whether there is evidence of genetic overlap (i.e., shared genes) more than expected by chance.'	
	dat4 = rbindlist(lapply(res, '[[', 4), fill=TRUE)
	dat4 = cbind(dat4[,c(1,2,4)], as.data.table(apply(dat4[,-c(1,2,4)], 2, as.numeric)))	
	format_header(wb, "T4.sptwas_cofdr", c(1, 14), h4, header_style)
	format_header(wb, "T4.sptwas_cofdr", c(2, 14), subh4, subheader_style, heights=3*13)
	format_main(wb, "T4.sptwas_cofdr", dat4, title_style, main_style)
	
	h5 = "Table 5. PWAS-cofdr analysis results for all validated genes across pairwise colocalization analyses"
	subh5 = 'Comparisons=trait comparisons identifying the "Gene"; Gene=HGNC symbol of the colocalized gene; TT=posterior possibility that corresponding Gene is shared between trait1 and trait2; Z_trait1=z-scores of the Gene estimated by PWAS for trait1; Z_trait2=z-scores of the Gene estimated by PWAS using FUSION for trait2; prop=the overall proportion of genes that are shared across both traits; LRT=the p-value of the likelihood-ratio test indicating whether there is evidence of genetic overlap (i.e., shared genes) more than expected by chance.'	
	dat5 = rbindlist(lapply(res, '[[', 5), fill=TRUE)
	dat5 = cbind(dat5[,1:2], as.data.table(apply(dat5[,-c(1:2)], 2, as.numeric)))	
	format_header(wb, "T5.pwas_cofdr", c(1, 14), h5, header_style)
	format_header(wb, "T5.pwas_cofdr", c(2, 14), subh5, subheader_style, heights=2*13)
	format_main(wb, "T5.pwas_cofdr", dat5, title_style, main_style)
	
	saveWorkbook(wb, file = save_path, overwrite = TRUE)
	
	return(save_path)
}


# Cofdr_Excel_Summary(excel_dir=paste0("03_cofdr/10_Excel_Summary/Tables/", c("A2", "B2", "C2"), "/respiratory"), save_path="10_Excel_Summary/Tables/A2_B2_C2_respiratory_more_info.xlsx")

# # generate Table 2.4 in COVID manuscript
# res = read_all_sheets("/mnt/data/xue/Data/02_COVID/03_cofdr/10_Excel_Summary/Tables/A2_B2_C2_respiratory_more_info.xlsx")
# dat = res$T4.spwas_cofdr[-1, ] %>% row_to_names(row_number = 1)
# final = data.frame()
# for(t1 in covidhgi){
	# for(t2 in respiratory){
		# pheno = paste0(t1, "_", t2)
		# tmp = dat[dat$Comparisons==pheno,]
		# tmp = tmp[tmp$FDR_trait1<0.05&tmp$FDR_trait2>0.05, ]
		# final = rbind(final, tmp[tmp$tissue=="Lung",], tmp[tmp$tissue=="Whole_Blood",], tmp[tmp$tissue=="smultixcan",])	
	# }
# }
# fwrite(final, file="tmp.txt", sep="\t")
