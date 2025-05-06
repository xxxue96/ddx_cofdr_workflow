library(openxlsx)
library(data.table)
source("03_cofdr/Cofdr_Workflow_params.r")
source("Excel_Format_Params.r")

# Apply this function if there is no sig.cofdr.snp 
Cofdr_Format_Excel_nosig <- function(t1, t2, save_dir)
{
	pheno = paste0(t1, "_", t2); print(paste(pheno, "start!"))
	gene2func = paste0("03_cofdr/05_02_Cofdr_FUMA_Result/GENE2FUNC/",pheno)
	unzip(paste0(gene2func,".zip"), exdir=gene2func)
	
	# create workbook
	wb <- createWorkbook()
	addWorksheet(wb, "T8.gene_based_cofdr")
	addWorksheet(wb, "T9.sig_gene_GSEA")
	addWorksheet(wb, "F1.sig_gene_DEG")

	# T8.gene_based_cofdr
	if(file.exists(paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno, "/", pheno, "_geneCompare.txt.gz"))){
		h8 = "Table 8. Gene-based co-fdr analysis to detect genes associated with both traits using MAGMA, PrediXcan or FUSION"
		subh8 = "Significantly colocalized genes were identified by integrating gene genomic positions, mRNA expression, splicing and plasma protein data with genome-wide summary statistics of COVID-19. Gene-based z-scores estimated through MAGMA, TWAS, spTWAS, and PWAS approaches were listed in the table, which were served as cofdr inputs to calculate corresponding posterior possibility (TT). Genes identified by at least three approaches were highlighted (total_detect>=3)."
		dat8 = fread(paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno, "/", pheno, "_geneCompare.txt.gz"))
		dat8 = dat8[order(dat8$total_detect, decreasing=TRUE), ] 
		format_header(wb, "T8.gene_based_cofdr", c(1, ncol(dat8)), h8, header_style)
		format_header(wb, "T8.gene_based_cofdr", c(2, ncol(dat8)), subh8, subheader_style, heights=13*2)
		format_main(wb, "T8.gene_based_cofdr", dat8, title_style, main_style)
		addStyle(wb, "T8.gene_based_cofdr", other_style, rows = which(dat8$total_detect>=3)+3, cols = 1, gridExpand = TRUE)
	} else{
		removeWorksheet(wb, "T8.gene_based_cofdr")
	}
	
	# T9.sig_gene_GSEA
	if(file.exists(paste0(gene2func,"/GS.txt"))){
		h9 = "Table 9. Enrichment of prioritized genes in Gene Sets"
		subh9 = "For colocalized genes highlighted in Table 8(total_detect>=3), significant gene sets with adjusted P-value<0.05 were shown here. Category=One of the category from MsigDB; GeneSet=Name of gene set as provided by MsigDB; N_genes=Number of genes in a gene set; N_overlap=Number of input genes overlapping with the gene set; p=Hypergeometric test (upper tail) P-value; adjP=Benjamini-Hochberg adjusted P-value; genes=Genes overlapping with the gene set; link=Link to the MsigDB page if available"
		dat9 = fread(paste0(gene2func,"/GS.txt"))
		dat9 = dat9[order(dat9$adjP, decreasing=FALSE),]
		format_header(wb, "T9.sig_gene_GSEA", c(1, 11), h9, header_style)
		format_header(wb, "T9.sig_gene_GSEA", c(2, 11), subh9, subheader_style, heights=13*4)
		format_main(wb, "T9.sig_gene_GSEA", dat9, title_style, main_style)
	} else{
		removeWorksheet(wb, "T9.sig_gene_GSEA")
	}
	
	# F1.sig_gene_DEG
	if(file.exists(paste0(gene2func,".png"))){
		figure1 = paste0(gene2func,".png")
		fh1 = "F1. Enrichment of prioritized genes in DEG Sets"
		subfh1 = "For colocalized genes highlighted in Table 8 (total_detect>=3), significantly enriched DEG sets (Pbon < 0.05) are highlighted in red."
		format_header(wb, "F1.sig_gene_DEG", c(1, 10), fh1, header_style)
		format_header(wb, "F1.sig_gene_DEG", c(2, 10), subfh1, subheader_style)
		insertImage(wb, "F1.sig_gene_DEG", figure1, startRow = 4, startCol = 1, width = 6.21, height = 5)
	} else{
		
		removeWorksheet(wb, "F1.sig_gene_DEG")
	}
			
	if(!identical(names(wb), character(0)))
	{
		save_path = paste0(save_dir, "/", pheno, ".xlsx")
		saveWorkbook(wb, file = save_path, overwrite = TRUE)
	}

	unlink(gene2func, recursive = TRUE)
	print(paste(pheno, "done!"))
}

# Apply this function if there is at least 1 sig.cofdr.snp
Cofdr_Format_Excel_sig <- function(t1, t2, save_dir, multi_cofdr=FALSE)
{
	pheno = paste(c(t1, t2), collapse="_"); print(paste(pheno, "start!"))
	header = as.vector(outer(c("FRQ_", "BETA_", "SE_", "P_", "Z_"), c(t1, t2), paste0))
	snp2gene = paste0("03_cofdr/05_02_Cofdr_FUMA_Result/SNP2GENE/",pheno)
	unzip(paste0(snp2gene,".zip"), exdir=snp2gene)
	gene2func = paste0("03_cofdr/05_02_Cofdr_FUMA_Result/GENE2FUNC/",pheno)
	unzip(paste0(gene2func,".zip"), exdir=gene2func)
	fuma_scale = 1e-5/(1-TTthres)
	
	# create workbook
	wb <- createWorkbook()
	addWorksheet(wb, "T1.sig_snps")
	addWorksheet(wb, "T2.risk_loci")
	addWorksheet(wb, "T3.GWAS_catalog")
	addWorksheet(wb, "T4.candidate_snps")
	addWorksheet(wb, "T5.annov_enrichment")
	addWorksheet(wb, "T6.mapping_genes")
	addWorksheet(wb, "T7.commonLoci_cofdr_gwaspw_hypr")
	addWorksheet(wb, "T8.gene_based_cofdr")
	addWorksheet(wb, "T9.sig_gene_GSEA")
	addWorksheet(wb, "F1.sig_gene_DEG")

	# T1.sig_snps
	h1 = paste("Table 1. Significantly colocalized SNPs detected by co-fdr at posterior possibility threshold", TTthres)
	
	if(multi_cofdr==TRUE){
		subh1 = "TT=posterior possibility that corresponding SNP is shared among multiple traits;A1=non-effect allele; A2=effect allele"
	} else{
		subh1 = "TT=posterior possibility that corresponding SNP is shared between trait1 and trait2;A1=non-effect allele; A2=effect allele"
	}
	
	dat1 = fread(paste0("03_cofdr/02_Cofdr_SigSNP/", pheno, "_", TTthres, "TT.txt.gz"))	
	dat1 = dat1[, c("SNP", "CHR", "BP", "A1", "A2", "TT", header), with=FALSE]
	dat1 = dat1[order(dat1$TT, decreasing=TRUE), ]
	format_header(wb, "T1.sig_snps", c(1, ncol(dat1)), h1, header_style)
	format_header(wb, "T1.sig_snps", c(2, ncol(dat1)), subh1, subheader_style)
	format_main(wb, "T1.sig_snps", dat1, title_style, main_style)

	# T2.risk_loci
	h2 = "Table 2. Identification of genomic risk loci"
	subh2 = "Genomic locus=Index of genomic risk loci; rsID=lead SNP with the largest TT in the genomic locus; nSNPS=number of candidate SNPs from input GWAS or reference panel in the genomic locus; nGWASSNPS=number of candidate SNPs from input GWAS; nIndSigSNPs=number of independent significant SNPs; IndSigSNPs=rsID of the independent significant SNPs in the genomic locus at r2 threshold 0.05; nLeadSNPs=number of lead SNPs; LeadSNPs=rsID of lead SNPs in the genomic locus at r2 threshold 0.05"
	dat2 = fread(paste0(snp2gene, "/GenomicRiskLoci.txt"))	
	dat2$TT = 1-dat2$p/fuma_scale
	dat2 = dat2[, -c("uniqID", "p")]
	dat2 = dat2[order(dat2$GenomicLocus, decreasing=FALSE), c(1:4,13,5:12)]
	format_header(wb, "T2.risk_loci", c(1, ncol(dat2)), h2, header_style)
	format_header(wb, "T2.risk_loci", c(2, ncol(dat2)), subh2, subheader_style, heights=13*3)
	format_main(wb, "T2.risk_loci", dat2, title_style, main_style)

	# T3.GWAS_catalog
	h3 = "Table 3. GWAS catalog annotations for SNPs in LD with independent significant SNPs in the genomic risk loci"
	subh3 = 'GenomicLocus=Index of genomic risk loci; IndSigSNP=rsID of the independent significant SNP in LD with the "snp" reported in GWAS catalog; TT_IndSigSNP=posterior possibility that the corresponding "IndSigSNP" is colocalized; chr=chromosome of IndSigSNP; bp=position on hg19 of IndSigSNP; snp=rsID of reported SNP in GWAS catalog; TT_snp=posterior possibility that the corresponding "snp" is colocalized; Trait=The trait reported in GWAScatalog; PMID=PubMed ID; Study=Title of paper; InitialN=Sample size and ancestry information for reported GWAS in GWAS catalog; ReplicationN=Sample size and ancestry information for replicated GWAS in GWAS catalog; Region=Cytogenetic region for the "snp"; ReportedGene=Associated genes of the "snp" reported by author; MappedGene=Summary of mapped genes or nearest upstream and downstream genes for the "snp"; UpGene=Nearest upstream genes if the "snp" not within the gene; DownGene=Nearest downstream genes if the "snp" not within the gene; SNP_Gene_ID=Mapped genes if the "snp" within the gene; UpGeneDist=distance in kb for "UpGene" to the "snp", if not within gene; DownGeneDist=distance in kb for "DownGene" to the "snp", if not within gene; Strongest=The most strongly associated "snp" with the trait and corresponding risk allele; SNPs="Strongest" SNPs; merged=Whether the "Strongest" SNP is merged to a subsequent rs record; SNP_ID_cur=Current rs number of the "Strongest" SNP; Content=The most severe consequence of the "Strongest" SNP on the "MappedGene"; Intergenic=Whether the "Strongest" SNP is intergenic or not; RiskAF=Risk allele frequency for the "Strongest" SNP; P=Reported p-value for the "Strongest" SNP; Pmlog=-log(P); Ptext=Information describing context of "P"; OrBeta=Reported OR or beta of the "Strongest" SNP risk allele; 95CI=Reported 95% confidence interval of the "Strongest" SNP risk allele; Platform=Genotyping platform for reported GWAS in GWAS catalog; CNV=Whether copy number variation study was performed or not'
	source("/mnt/data/xue/Data/02_COVID/04_script/08_Annot_gwas_catalog.r")
	annot_path = paste0("03_cofdr/01_Cofdr_Result/", pheno, ".txt.gz")
	catalog_path = paste0(snp2gene, "/gwascatalog.txt")
	Annot_gwas_catalog(annot_path, catalog_path, save_path=catalog_path)
	dat3 = fread(catalog_path)
	dat3 = dat3[order(dat3$TT_snp, decreasing=TRUE),]
	format_header(wb, "T3.GWAS_catalog", c(1, ncol(dat3)), h3, header_style)
	format_header(wb, "T3.GWAS_catalog", c(2, ncol(dat3)), subh3, subheader_style, heights=13*4)
	format_main(wb, "T3.GWAS_catalog", dat3, title_style, main_style)

	# T4.candidate_snps
	h4 = "Table 4. Functional annotations of candidate SNPs in LD of independent significant SNPs at r2 theshold 0.05"
	subh4 = 'Candidate SNPs were prioritized by CADD score and those with CADD scores not less than 12.37 were highlighted. rsID=rsID for annotated candidate SNPs; TT=posterior possibility of corresponding SNPs from GWAS. Summary statistics of SNPs extracted from the reference panel are denoted as "NA"; r2=r2 of the candidate SNP in "rsID" with the independent significant SNP in "IndSigSNP"; Genomic locus=Index of genomic risk loci shown in table "T2.risk_loci"; nearestGene=The nearest gene of the SNP; dist=Distance between the SNP and its nearest gene; func=Functional consequence of the SNP on the gene; CADD=CADD score indicating SNP deleteriousness. SNPs with higher CADD scores tend to be more deleterious; RDB=RegulomeDB categorical score. The higher the score (form 1a to 7), the more evidence for the SNP having regulatory functions; minChrState/commonChrState=The minimum/most common chromatin state, a marker of chromatin accessibility. The lower the score, the more open the chromatin is. posMapFilt/eqtlMapFilt/ciMapFilt=Whether the SNP was included in positional mapping/eQTL mapping/chromatin interaction mapping or not. 1 is used and 0 is not.'
	dat4 = fread(paste0(snp2gene, "/snps.txt"))
	dat4$TT = 1-dat4$gwasP/fuma_scale
	dat4 = dat4[, -c("uniqID", "gwasP")]	
	dat4 = dat4[order(dat4$CADD, decreasing=TRUE), c(1:6,20,7:19)]
	format_header(wb, "T4.candidate_snps", c(1, ncol(dat4)), h4, header_style)
	format_header(wb, "T4.candidate_snps", c(2, ncol(dat4)), subh4, subheader_style, heights=13*5)
	format_main(wb, "T4.candidate_snps", dat4, title_style, main_style)
	addStyle(wb, "T4.candidate_snps", other_style, rows = which(dat4$CADD>12.37)+3, cols = 1, gridExpand = TRUE)

	# T5.annov_enrichment
	h5 = "Table 5. Enrichment test of functional consequences of candidate SNPs on genes"
	subh5 = "ANNOVAR enrichment test of functional consequences of candidate SNPs against the 1000G reference panel, and significant results are highlighted(fisher.P<0.05). Consequences=Functional effects of SNPs; ref.count=Number of SNPs in the reference panel annotated with the corresponding functional category; ref.prop=The proportion of annotated SNPs within all SNPs in the reference panel; count=Number of candidate SNPs annotated with the corresponding functional category; prop=The proportion of annotated candidate SNPs; enrichment=prop/ref.prop, a value to show whether the annotation is enriched (larger than 1) or depleted otherwise; fisher.P=two-sided p-value from Fisher’s exact test, a statistic to show whether enrichment of functional consequences of candidate SNPs is significantly different from the reference panel or not."
	dat5 = fread(paste0(snp2gene, "/annov.stats.txt"))
	dat5 = dat5[order(dat5$fisher.P, decreasing=FALSE), ]
	format_header(wb, "T5.annov_enrichment", c(1, 14), h5, header_style)
	format_header(wb, "T5.annov_enrichment", c(2, 14), subh5, subheader_style, heights=13*5)
	format_main(wb, "T5.annov_enrichment", dat5, title_style, main_style)
	addStyle(wb, "T5.annov_enrichment", other_style, rows = which(dat5$fisher.P<0.05)+3, cols = 1, gridExpand = TRUE)

	# T6.mapping_genes
	h6 = "Table 6. Gene mapping results based on positional, eQTL and chromatin interaction"
	subh6 = 'Genes with SNPs mapped by positional, eQTL, and chromatin interaction were highlighted. pL1=The probability of a gene to be loss-of-function intolerant, and the higher score indicates the gene is more intolerant; ncRVIS=The score of a gene to be intolerant to non-coding residual variations, and the higher score indicates the gene is more intolerant; posMap/eqtlMap/ciMap=Whether there are any SNPs mapped to the gene by positional mapping/eQTL mapping/chromatin interaction mapping or not; posMapSNPs/eqtlMapSNPs=Number of mapped SNPs by positional mapping/eQTL mapping; posMapMaxCADD=The maximum CADD score of mapped SNPs by positional mapping; eqtlMapminP/eqtlMapminQ=The minimum p-value/FDR for mapped SNPs by eQTL mapping; eqtlMapts/ciMapts=Tissue types to detect significant eQTL associations/chromatin interaction between mapped SNPs and the gene; eqtlDirection=Direction of the effect of tested allele on gene expression, which can be "+" (tested allele increases gene expression) or "-" (tested allele decreases gene expression); minGwasP=The minimum p-value of SNPs mapped to the gene; IndSigSNPs=rsID of independent significant SNPs in LD with mapped SNPs; GenomicLocus=Index of genomic risk loci containing mapped SNPs.'
	dat6 = fread(paste0(snp2gene, "/genes.txt"))
	dat6 <- transform(dat6, posMap=ifelse(posMapSNPs==0, "No", "Yes"), eqtlMap=ifelse(eqtlMapSNPs==0, "No", "Yes"))
	dat6 = dat6[, c(1, 8:9, 3:7, 24:25, 19, 10:18, 20:23)]
	dat6 = dat6[order(dat6$posMap, dat6$eqtlMap, dat6$ciMap, decreasing=TRUE), ]
	format_header(wb, "T6.mapping_genes", c(1, ncol(dat6)), h6, header_style)
	format_header(wb, "T6.mapping_genes", c(2, ncol(dat6)), subh6, subheader_style, heights=13*4)
	format_main(wb, "T6.mapping_genes", dat6, title_style, main_style)
	addStyle(wb, "T6.mapping_genes", other_style, rows = which(dat6$posMap=="Yes" & dat6$eqtlMap=="Yes" & dat6$ciMap=="Yes")+3, cols = 1, gridExpand = TRUE)

	# T7.commonLoci_cofdr_gwaspw_hypr
	if(multi_cofdr==TRUE){
		h7 = "Table 7.Significantly colocalized SNPs detected by cofdr and hyprcoloc"
		subh7 = 'Significant SNPs detected by cofdr and hyprcoloc were highlighted (total_detect=2). TT_cofdr=posterior possibility from co-fdr that the corresponding SNP in the "SNP" column is colocalized; genomic region is the colocalized SNP; posterior_prob_hyprcoloc=posteior possibility from hyprcoloc that a defined LD-independent region contains a colocalized SNP; posterior_explained_by_snp_hyprcoloc=posterior possibility from hyprcoloc that the corresponding SNP located in the defined genomic region is the colocalized SNP.'
	} else{
		h7 = "Table 7.Significantly colocalized SNPs detected by cofdr, gwas-pw and hyprcoloc"
		subh7 = 'Significant SNPs detected by cofdr, gwas-pw and hyprcoloc were highlighted (total_detect=3). TT_cofdr=posterior possibility from co-fdr that the corresponding SNP in the "SNP" column is colocalized; chunkPPA_3_gwaspw=posteior possibility from gwas-pw that a defined LD-independent region contains a colocalized SNP; snpPPA_3_gwaspw=posterior possibility from gwas-pw that the corresponding SNP located in the defined genomic region is the colocalized SNP; posterior_prob_hyprcoloc=posteior possibility from hyprcoloc that a defined LD-independent region contains a colocalized SNP; posterior_explained_by_snp_hyprcoloc=posterior possibility from hyprcoloc that the corresponding SNP located in the defined genomic region is the colocalized SNP.'
	}
	
	if(file.exists(paste0("03_cofdr/09_sigSNP_Comparison/", pheno, ".txt.gz"))){		
		dat7 = fread(paste0("03_cofdr/09_sigSNP_Comparison/", pheno, ".txt.gz"))
		colnames(dat7)[c(2,4,5,6,7)] = c("TT_cofdr", "chunkPPA_3_gwaspw", "snpPPA_3_gwaspw", "posterior_prob_hyprcoloc", "posterior_explained_by_snp_hyprcoloc")
		dat7 = dat7[, c("SNP", "CHR", "BP", "A1", "A2", header, "TT_cofdr", "chunkPPA_3_gwaspw", "snpPPA_3_gwaspw", "posterior_prob_hyprcoloc", "posterior_explained_by_snp_hyprcoloc", "cofdr_detect", "gwaspw_detect", "hyprcoloc_detect", "total_detect"), with=FALSE]		
		dat7 = dat7[order(dat7$total_detect, decreasing=TRUE), ] 
	} else{
		tmp = c("SNP", "CHR", "BP", "A1", "A2", header, "TT_cofdr", "chunkPPA_3_gwaspw", "snpPPA_3_gwaspw", "posterior_prob_hyprcoloc", "posterior_explained_by_snp_hyprcoloc", "cofdr_detect", "gwaspw_detect", "hyprcoloc_detect", "total_detect")
		dat7 = setNames(data.frame(matrix(ncol = length(tmp), nrow = 0)), tmp)		
	}	
	format_header(wb, "T7.commonLoci_cofdr_gwaspw_hypr", c(1, ncol(dat7)), h7, header_style)
	format_header(wb, "T7.commonLoci_cofdr_gwaspw_hypr", c(2, ncol(dat7)), subh7, subheader_style, heights=13*3)
	format_main(wb, "T7.commonLoci_cofdr_gwaspw_hypr", dat7, title_style, main_style)
	addStyle(wb, "T7.commonLoci_cofdr_gwaspw_hypr", other_style, rows = which(dat7$total_detect==3)+3, cols = 1, gridExpand = TRUE)
	if(multi_cofdr==TRUE){renameWorksheet(wb, "T7.commonLoci_cofdr_gwaspw_hypr", "T7.commonLoci_cofdr_hypr")}

	# T8.gene_based_cofdr
	if(file.exists(paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno, "/", pheno, "_geneCompare.txt.gz"))){
		h8 = "Table 8. Gene-based co-fdr analysis to detect genes associated with both traits using MAGMA, PrediXcan or FUSION"
		subh8 = "Significantly colocalized genes were identified by integrating gene genomic positions, mRNA expression, splicing and plasma protein data with genome-wide summary statistics of COVID-19. Gene-based z-scores estimated through MAGMA, TWAS, spTWAS, and PWAS approaches were listed in the table, which were served as cofdr inputs to calculate corresponding posterior possibility (TT). Genes identified by at least three approaches were highlighted (total_detect>=3)."
		dat8 = fread(paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno, "/", pheno, "_geneCompare.txt.gz"), fill=TRUE)
		dat8 = dat8[order(dat8$total_detect, decreasing=TRUE), ] 
		format_header(wb, "T8.gene_based_cofdr", c(1, ncol(dat8)), h8, header_style)
		format_header(wb, "T8.gene_based_cofdr", c(2, ncol(dat8)), subh8, subheader_style, heights=13*2)
		format_main(wb, "T8.gene_based_cofdr", dat8, title_style, main_style)
		addStyle(wb, "T8.gene_based_cofdr", other_style, rows = which(dat8$total_detect>=3)+3, cols = 1, gridExpand = TRUE)
	} else{
		removeWorksheet(wb, "T8.gene_based_cofdr")
	}
	
	# T9.sig_gene_GSEA
	if(file.exists(paste0(gene2func,"/GS.txt"))){
		h9 = "Table 9. Enrichment of prioritized genes in Gene Sets"
		subh9 = "For colocalized genes highlighted in Table 8(total_detect>=3), significant gene sets with adjusted P-value<0.05 were shown here. Category=One of the category from MsigDB; GeneSet=Name of gene set as provided by MsigDB; N_genes=Number of genes in a gene set; N_overlap=Number of input genes overlapping with the gene set; p=Hypergeometric test (upper tail) P-value; adjP=Benjamini-Hochberg adjusted P-value; genes=Genes overlapping with the gene set; link=Link to the MsigDB page if available"
		dat9 = fread(paste0(gene2func,"/GS.txt"))
		dat9 = dat9[order(dat9$adjP, decreasing=FALSE),]
		format_header(wb, "T9.sig_gene_GSEA", c(1, 11), h9, header_style)
		format_header(wb, "T9.sig_gene_GSEA", c(2, 11), subh9, subheader_style, heights=13*4)
		format_main(wb, "T9.sig_gene_GSEA", dat9, title_style, main_style)
	} else{
		removeWorksheet(wb, "T9.sig_gene_GSEA")
	}
	
	# F1.sig_gene_DEG
	if(file.exists(paste0(gene2func,".png"))){
		figure1 = paste0(gene2func,".png")
		fh1 = "F1. Enrichment of prioritized genes in DEG Sets"
		subfh1 = "For colocalized genes highlighted in Table 8 (total_detect>=3), significantly enriched DEG sets (Pbon < 0.05) are highlighted in red."
		format_header(wb, "F1.sig_gene_DEG", c(1, 10), fh1, header_style)
		format_header(wb, "F1.sig_gene_DEG", c(2, 10), subfh1, subheader_style)
		insertImage(wb, "F1.sig_gene_DEG", figure1, startRow = 4, startCol = 1, width = 6.21, height = 5)
	} else{		
		removeWorksheet(wb, "F1.sig_gene_DEG")
	}
	
	save_path = paste0(save_dir, "/", pheno, ".xlsx")
	saveWorkbook(wb, file = save_path, overwrite = TRUE)

	unlink(snp2gene, recursive = TRUE)
	unlink(gene2func, recursive = TRUE)
	print(paste(pheno, "done!"))
}

Cofdr_Format_Excel <- function(t1, t2, save_dir, multi_cofdr=FALSE)
{
	pheno = paste(c(t1, t2), collapse="_")

	if(file.exists(paste0("03_cofdr/05_02_Cofdr_FUMA_Result/SNP2GENE/", pheno, ".zip"))){
		Cofdr_Format_Excel_sig(t1, t2, save_dir, multi_cofdr)
	} else{
		Cofdr_Format_Excel_nosig(t1, t2, save_dir)
	}
	
}

# covidhgi = c("A2", "B2", "C2")
# respiratory = c("Asthma", "ILD", "IPF", "COPD", "Pneumonia.meta", "Pneumonia.b.meta", "Pneumonia.v.meta", "Flu.meta")

# for(t1 in covidhgi)
# {
	# for(t2 in respiratory){ Cofdr_Format_Excel(t1, t2, save_dir = paste0("03_cofdr/10_Excel_Summary/Tables/", t1, "/respiratory")) }
	
# }

# Cofdr_Format_Excel(t1="A2", t2=c("Asthma", "ILD", "IPF", "COPD"), save_dir="11_Multi_Cofdr_Summary/Tables/A2", multi_cofdr=TRUE)
# Cofdr_Format_Excel(t1="B2", t2=c("Asthma", "ILD", "IPF"), save_dir="11_Multi_Cofdr_Summary/Tables/B2", multi_cofdr=TRUE)
# Cofdr_Format_Excel(t1="C2", t2=c("Asthma", "ILD", "IPF"), save_dir="11_Multi_Cofdr_Summary/Tables/C2", multi_cofdr=TRUE)