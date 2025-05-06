library(openxlsx)
library(data.table)
library(dplyr)
source("02_ddx/DDx_Workflow_params.r")
source("Excel_Format_Params.r")

DDx_Format_Excel <- function(t1, t2, save_dir)
{
	pheno = paste0(t1, "_", t2); print(paste(pheno, "start!"))
	snp2gene = paste0("02_ddx/04_FUMA_Results/SNP2GENE/",pheno)
	snp2gene_t1 = paste0("01_data/04_FUMA_Results/SNP2GENE/",t1)
	snp2gene_t2 = paste0("01_data/04_FUMA_Results/SNP2GENE/",t2)
	unzip(paste0(snp2gene,".zip"), exdir=snp2gene)
	unzip(paste0(snp2gene_t1,".zip"), exdir=snp2gene_t1)
	unzip(paste0(snp2gene_t2,".zip"), exdir=snp2gene_t2)
	celltype = paste0("02_ddx/04_FUMA_Results/Cell_Type/",pheno)
	celltype_t1 = paste0("01_data/04_FUMA_Results/Cell_Type/",t1)
	celltype_t2 = paste0("01_data/04_FUMA_Results/Cell_Type/",t2)
	unzip(paste0(celltype,".zip"), exdir=celltype)
	unzip(paste0(celltype_t1,".zip"), exdir=celltype_t1)
	unzip(paste0(celltype_t2,".zip"), exdir=celltype_t2)
	image_dir = "02_ddx/04_FUMA_Results/Image/05_Cell_Type"

	# create workbook
	wb <- createWorkbook()
	addWorksheet(wb, "T1.sig_snps")
	addWorksheet(wb, "T2.risk_loci")
	addWorksheet(wb, "T3.GWAS_catalog")
	addWorksheet(wb, "T4.candidate_snps")
	addWorksheet(wb, "T5.annov_enrichment")
	addWorksheet(wb, "T6.mapping_genes")
	addWorksheet(wb, "T7.commonLoci_ddx_ccgwas")
	addWorksheet(wb, "T8.gene_ddx")
	addWorksheet(wb, "T9.gene_ccgwas")
	addWorksheet(wb, "T10.commonGene_ddx_ccgwas")
	addWorksheet(wb, "T11.magma_gene_set")
	addWorksheet(wb, "T12.magma_tissue")
	addWorksheet(wb, "T13.magma_cell_type")
	addWorksheet(wb, "F1.magma_cell_type")
		
	# T1.sig_snps
	h1 = paste("Table 1. Significantly differentiated SNPs detected by DDx at p-value threshold", pthres)
	subh1 = "A1: non-effect allele; A2: effect allele"
	dat1 = fread(paste0("02_ddx/02_DDx_SigSNP/", pheno, "_", pthres, ".txt.gz"))
	dat1 = dat1[, colnames(dat1)[!grepl("^N_", colnames(dat1))], with=FALSE]
	dat1 = dat1[order(dat1$P, decreasing=FALSE), ]
	format_header(wb, "T1.sig_snps", c(1, ncol(dat1)), h1, header_style)
	format_header(wb, "T1.sig_snps", c(2, ncol(dat1)), subh1, subheader_style)
	format_main(wb, "T1.sig_snps", dat1, title_style, main_style)

	# T2.risk_loci
	h2 = "Table 2. Identification of genomic risk loci"
	subh2 = "Genomic locus=Index of genomic risk loci; rsID=lead SNP with the smalllest p-value in the genomic locus; nSNPS=number of candidate SNPs from input GWAS or reference panel in the genomic locus; nGWASSNPS=number of candidate SNPs from input GWAS; nIndSigSNPs=number of independent significant SNPs; IndSigSNPs=rsID of the independent significant SNPs in the genomic locus at r2 threshold 0.05; nLeadSNPs=number of lead SNPs; LeadSNPs=rsID of lead SNPs in the genomic locus at r2 threshold 0.05"
	dat2 = fread(paste0(snp2gene, "/GenomicRiskLoci.txt"))
	dat2 = dat2[, -"uniqID"]
	dat2 = dat2[order(dat2$GenomicLocus, decreasing=FALSE), ]
	format_header(wb, "T2.risk_loci", c(1, ncol(dat2)), h2, header_style)
	format_header(wb, "T2.risk_loci", c(2, ncol(dat2)), subh2, subheader_style, heights=13*4)
	format_main(wb, "T2.risk_loci", dat2, title_style, main_style)

	# T3.GWAS_catalog
	h3 = "Table 3. GWAS catalog annotations for SNPs in LD with independent significant SNPs in the genomic risk loci"
	subh3 = 'GenomicLocus=Index of genomic risk loci; IndSigSNP=rsID of the independent significant SNP in LD with the "snp" reported in GWAS catalog; P_IndSigSNP=posterior possibility that the corresponding "IndSigSNP" is colocalized; chr=chromosome of IndSigSNP; bp=position on hg19 of IndSigSNP; snp=rsID of reported SNP in GWAS catalog; P_snp=posterior possibility that the corresponding "snp" is colocalized; Trait=The trait reported in GWAScatalog; PMID=PubMed ID; Study=Title of paper; InitialN=Sample size and ancestry information for reported GWAS in GWAS catalog; ReplicationN= Sample size and ancestry information for replicated GWAS in GWAS catalog; Region=Cytogenetic region for the "snp"; ReportedGene=Associated genes of the "snp" reported by author; MappedGene=Summary of mapped genes or nearest upstream and downstream genes for the "snp"; UpGene=Nearest upstream genes if the "snp" not within the gene; DownGene=Nearest downstream genes if the "snp" not within the gene; SNP_Gene_ID=Mapped genes if the "snp" within the gene; UpGeneDist=distance in kb for "UpGene" to the "snp", if not within gene; DownGeneDist=distance in kb for "DownGene" to the "snp", if not within gene; Strongest=The most strongly associated "snp" with the trait and corresponding risk allele; SNPs="Strongest" SNPs; merged=Whether the "Strongest" SNP is merged to a subsequent rs record; SNP_ID_cur=Current rs number of the "Strongest" SNP; Content=The most severe consequence of the "Strongest" SNP on the "MappedGene"; Intergenic=Whether the "Strongest" SNP is intergenic or not; RiskAF=Risk allele frequency for the "Strongest" SNP; P=Reported p-value for the "Strongest" SNP; Pmlog=-log(P); Ptext=Information describing context of "P"; OrBeta=Reported OR or beta of the "Strongest" SNP risk allele; 95CI=Reported 95% confidence interval of the "Strongest" SNP risk allele; Platform=Genotyping platform for reported GWAS in GWAS catalog; CNV=Whether copy number variation study was performed or not'
	source("04_script/08_Annot_gwas_catalog.r")
	annot_path = paste0("02_ddx/01_DDx_Results/", pheno, ".txt.gz")
	catalog_path = paste0(snp2gene, "/gwascatalog.txt")
	Annot_gwas_catalog(annot_path, catalog_path, save_path=catalog_path)
	dat3 = fread(catalog_path)
	colnames(dat3)[c(3,7)] = c("P_IndSigSNP", "P_snp")
	dat3 = dat3[order(dat3$P_snp, decreasing=FALSE),]
	format_header(wb, "T3.GWAS_catalog", c(1, ncol(dat3)), h3, header_style)
	format_header(wb, "T3.GWAS_catalog", c(2, ncol(dat3)), subh3, subheader_style, heights=13*4)
	format_main(wb, "T3.GWAS_catalog", dat3, title_style, main_style)

	# T4.candidate_snps
	h4 = "Table 4. Functional annotations of candidate SNPs in LD of independent significant SNPs at r2 theshold 0.05"
	subh4 = 'Candidate SNPs were prioritized by CADD score and those with CADD scores larger than 12.37 were highlighted. rsID=rsID for annotated candidate SNPs; gwasP/or/beta/se=p-value/odds ratio/effect size/standard error of SNPs from GWAS. Summary statistics of SNPs extracted from the reference panel are denoted as NA; r2=r2 of the candidate SNP in "rsID" with the independent significant SNP in "IndSigSNP"; Genomic locus=Index of genomic risk loci shown in Table 2; nearestGene=The nearest gene of the SNP; dist=Distance between the SNP and its nearest gene; func=Functional consequence of the SNP on the gene; CADD=CADD score indicating SNP deleteriousness. SNPs with higher CADD scores tend to be more deleterious; RDB=RegulomeDB categorical score. The higher the score (form 1a to 7), the more evidence for the SNP having regulatory functions; minChrState/commonChrState=The minimum/most common chromatin state, a marker of chromatin accessibility. The lower the score, the more open the chromatin is. posMapFilt/eqtlMapFilt/ciMapFilt=Whether the SNP was included in positional mapping/eQTL mapping/chromatin interaction mapping or not. 1 is used and 0 is not.'
	dat4 = fread(paste0(snp2gene, "/snps.txt"))
	dat4 = dat4[, -"uniqID"]
	dat4 = dat4[order(dat4$CADD, decreasing=TRUE), ]
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
	
	# T7.commonLoci_ddx_ccgwas
	h7 = "Table 7. Significantly differentiated loci detected by DDx and CCGWAS"
	subh7 = 'GenomicLocus=number of risk loci commonly detected by both DDx and CCGWAS; SNP=lead SNP of the "GenomicLocus"; A1=non-effect allele; A2=effect allele; FRQ=effect allele frequency; BETA=effect size of "SNP" estimated by DDx; SE=standard error of "SE"; P=pvalue of "SNP" estimated by DDx; Z=zscore of "SNP" estimated by DDx; N=total sample size involved in DDx analysis; OLS_beta/OLS_se/OLS_pval=case-case associations estimated by the CCGWAS OLS component; Exact_beta/Exact_se/Exact_pval=case-case associations estimated by the CCGWAS Exact component. The last 10 columns refer to the original GWAS statistics for two compared traits respectively.'
	dat7 = data.frame(matrix(ncol = 28, nrow = 0))		
	colnames(dat7) =  c("GenomicLocus", "SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "Z", "N", paste0("OLS_",c("beta","se","pval")), paste0("Exact_",c("beta","se","pval")), as.vector(outer(c("FRQ_", "BETA_", "SE_", "P_", "Z_"), c(t1,t2), paste0)))
	if(file.exists(paste0("02_ddx/07_SNP_Compare/", pheno, ".txt.gz"))){
		commonSNP = fread(paste0("02_ddx/07_SNP_Compare/", pheno, ".txt.gz"))[, colnames(dat7)[-1], with=FALSE]
		if(nrow(commonSNP)>0){		
			dat7 = merge(commonSNP, dat2[, "rsID"], by.x="SNP", by.y="rsID")
			dat7 = cbind(GenomicLocus=1:nrow(dat7), dat7)		
		}
	} 
	format_header(wb, "T7.commonLoci_ddx_ccgwas", c(1, ncol(dat7)), h7, header_style)
	format_header(wb, "T7.commonLoci_ddx_ccgwas", c(2, ncol(dat7)), subh7, subheader_style, heights=13*2)
	format_main(wb, "T7.commonLoci_ddx_ccgwas", dat7, title_style, main_style)
	
	# T8.gene_ddx
	h8 = "Table 8. Results summary of MAGMA/TWAS/spTWAS/PWAS-DDx"
	subh8 = 'Significant differentially associated genes were identified by integrating gene genomic positions, mRNA expression, splicing and plasma protein data with DDx-derived GWAS using MAGMA gene analysis, PrediXcan or FUSION. P_DDx_magma=the gene p-value estimated by MAGMA based on DDx-derived GWAS; FDR_DDx_magma=adjusted P_DDx_magma; P_DDx_twas=the gene p-value estimated by TWAS based on DDx-derived GWAS; tissue_twas=the eQTL model of every single tissue used to obtain gene-based z-scores in TWAS analysis based on DDx-derived GWAS; FDR_DDx_twas=adjusted P_DDx_twas; Z_DDx_twas=gene-based z-scores in TWAS analysis; P_DDx_sptwas=the gene p-value estimated by spTWAS based on DDx-derived GWAS; tissue_sptwas=the sQTL model of every single tissue used to obtain gene-based z-scores in spTWAS analysis based on DDx-derived GWAS; FDR_DDx_sptwas=adjusted P_DDx_sptwas; Z_DDx_sptwas=gene-based z-scores in spTWAS analysis; P_DDx_pwas=the gene p-value estimated by PWAS based on DDx-derived GWAS; FDR_DDx_pwas=adjusted P_DDx_pwas; Z_DDx_pwas=gene-based z-scores in PWAS analysis. Other columns refer to MAGMA/TWAS/spTWAS/PWAS estimates based on original GWAS of analyzed traits.'
	dat8 = fread(paste0("02_ddx/08_DDx_Gene_Result/", pheno, "/", pheno, "_geneCompare.txt.gz"))
	dat8 = dat8[order(dat8$total_detect, decreasing=TRUE), ]
	colnames(dat8)[1] = "Gene"
	format_header(wb, "T8.gene_ddx", c(1, ncol(dat8)), h8, header_style)
	format_header(wb, "T8.gene_ddx", c(2, ncol(dat8)), subh8, subheader_style, heights=13*3)
	format_main(wb, "T8.gene_ddx", dat8, title_style, main_style)
	
	# T9.gene_ccgwas
	h9 = "Table 9. Results summary of MAGMA/TWAS/spTWAS/PWAS-CC-GWAS"
	subh9 = 'Significant differentially associated genes were identified by integrating gene genomic positions, mRNA expression, splicing and plasma protein data withCC-GWAS OLS componnet using MAGMA gene analysis, PrediXcan or FUSION. P_OLS_magma=the gene p-value estimated by MAGMA based on CC-GWAS OLS componnet; FDR_OLS_magma=adjusted P_OLS_magma; P_OLS_twas=the gene p-value estimated by TWAS based on CC-GWAS OLS componnet; tissue_twas=the eQTL model of every single tissue used to obtain gene-based z-scores in TWAS analysis based on CC-GWAS OLS componnet; FDR_OLS_twas=adjusted P_OLS_twas; Z_OLS_twas=gene-based z-scores in TWAS analysis; P_OLS_sptwas=the gene p-value estimated by spTWAS based on CC-GWAS OLS componnet; tissue_sptwas=the sQTL model of every single tissue used to obtain gene-based z-scores in spTWAS analysis based on CC-GWAS OLS componnet; FDR_OLS_sptwas=adjusted P_OLS_sptwas; Z_OLS_sptwas=gene-based z-scores in spTWAS analysis; P_OLS_pwas=the gene p-value estimated by PWAS based on CC-GWAS OLS componnet; FDR_OLS_pwas=adjusted P_OLS_pwas; Z_OLS_pwas=gene-based z-scores in PWAS analysis. Other columns refer to MAGMA/TWAS/spTWAS/PWAS estimates based on original GWAS of analyzed traits.'
	if(file.exists(paste0("02_ddx/09_CCGWAS_Gene_Result/", pheno, "/OLS/", pheno, "_geneCompare.txt.gz"))){
		dat9 = fread(paste0("02_ddx/09_CCGWAS_Gene_Result/", pheno, "/OLS/", pheno, "_geneCompare.txt.gz"))
		dat9 = dat9[order(dat9$total_detect, decreasing=TRUE), ]
		colnames(dat9)[1] = "Gene"
	} else{
		dat9 = data.frame(matrix(nrow=0, ncol=44))
		colnames(dat9) = gsub("_DDx", "_OLS", colnames(dat8))
	}	
	format_header(wb, "T9.gene_ccgwas", c(1, ncol(dat9)), h9, header_style)
	format_header(wb, "T9.gene_ccgwas", c(2, ncol(dat9)), subh9, subheader_style, heights=13*3)
	format_main(wb, "T9.gene_ccgwas", dat9, title_style, main_style)

	# T10.commonGene_ddx_ccgwas
	h10 = "Table 10. Results summary of validated differential genes"
	subh10 = 'Results of MAGMA/TWAS/spTWAS/PWAS-DDx and MAGMA/TWAS/spTWAS/PWAS-CC-GWAS were merged respectively, then corresponding four groups of merged data were compared further and overlapping genes shown in at least three approaches were highlighted as validated differential genes (total_detect>=3). Gene-based p-value/FDR/z-scores estimated through MAGMA, TWAS, spTWAS, and PWAS approaches based on DDx/CC-GWAS OLS case-case associations and original GWAS data were listed in the table. magma_detect=1 indicates the "Gene" was detected by both MAGMA-DDx and MAGMA-CC-GWAS, and 0 otherwise; twas_detect=1 indicates the "Gene" was detected by both TWAS-DDx and TWAS-CC-GWAS, and 0 otherwise; sptwas_detect=1 indicates the "Gene" was detected by both spTWAS-DDx and spTWAS-CC-GWAS, and 0 otherwise; pwas_detect=1 indicates the "Gene" was detected by both PWAS-DDx and PWAS-CC-GWAS, and 0 otherwise; total_detect=magma_detect + twas_detect + sptwas_detect + pwas_detect'
	if(file.exists(paste0("02_ddx/10_DDx_CCGWAS_Gene_Compare/", pheno, "_geneCompare.txt.gz"))){
		dat10 = fread(paste0("02_ddx/10_DDx_CCGWAS_Gene_Compare/", pheno, "_geneCompare.txt.gz"))
		dat10 = dat10[order(dat10$total_detect, decreasing=TRUE), ]
		colnames(dat10)[1] = "Gene"
	} else{
		dat10 = data.frame(matrix(nrow=0, ncol=56))
		header = colnames(dat8)		
		header = append(header, paste0(c("P","FDR","Z"), "_OLS_magma"), grep("Z_DDx_magma", header))
		header = append(header, paste0(c("P","FDR","Z"), "_OLS_twas"), grep("Z_DDx_twas", header))
		header = append(header, paste0(c("P","FDR","Z"), "_OLS_sptwas"), grep("Z_DDx_sptwas", header))
		header = append(header, paste0(c("P","FDR","Z"), "_OLS_pwas"), grep("Z_DDx_pwas", header))		
		colnames(dat10) = header
	}
	format_header(wb, "T10.commonGene_ddx_ccgwas", c(1, ncol(dat10)), h10, header_style)
	format_header(wb, "T10.commonGene_ddx_ccgwas", c(2, ncol(dat10)), subh10, subheader_style, heights=13*2)
	format_main(wb, "T10.commonGene_ddx_ccgwas", dat10, title_style, main_style)
	addStyle(wb, "T10.commonGene_ddx_ccgwas", other_style, rows = which(dat10$total_detect>=3)+3, cols = 1, gridExpand = TRUE)
	
	# T11.magma_gene_set
	read_gene_set <- function(fuma_dir, trait)
	{		
		dat = fread(paste0("grep -v '^#' ", fuma_dir, "/magma.gsa.out"))	
		dat$FDR = p.adjust(dat$P, method = "fdr")
		dat = dat[, c("FULL_NAME", "P", "FDR")]
		colnames(dat) = c("Gene_set", paste0(c("P_", "FDR_"), trait))
		return(dat)
	}	
	h11 = "Table 11. Gene set analysis based on DDx-derived GWAS using MAGMA"
	subh11 = "Significantly associated gene sets at FDR threshold 0.05 were highlighted. FDR for each gene set was estimated by MAGMA based on DDx-derived GWAS or original GWAS of two compared traits"
	dat11 = read_gene_set(snp2gene, trait="DDx")
	dat11 = list(dat11, read_gene_set(snp2gene_t1, t1), read_gene_set(snp2gene_t2, t2)) %>% reduce(left_join)		
	dat11 = dat11[order(dat11$FDR_DDx, decreasing=FALSE), ]
	format_header(wb, "T11.magma_gene_set", c(1, ncol(dat11)), h11, header_style)
	format_header(wb, "T11.magma_gene_set", c(2, ncol(dat11)), subh11, subheader_style, heights=13*3)
	format_main(wb, "T11.magma_gene_set", dat11, title_style, main_style)
	addStyle(wb, "T11.magma_gene_set", other_style, rows = which(dat11$FDR_DDx<0.05)+3, cols = 1, gridExpand = TRUE)
	
	# T12.magma_tissue
	read_ts <- function(fuma_dir, trait)
	{
		dat = fread(paste0("grep -v '^#' ", fuma_dir, "/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out"))	
		dat$FDR = p.adjust(dat$P, method = "fdr")
		dat = dat[, c("FULL_NAME", "P", "FDR")]
		colnames(dat) = c("Tissue", paste0(c("P_", "FDR_"), trait))
		return(dat)
	}
	h12 = "Table 12. Tissue enrichment analysis based on DDx-derived GWAS using MAGMA"
	subh12 = "Significantly associated tissues at FDR threshold 0.05 were highlighted. FDR for each tissue was estimated by MAGMA based on DDx-derived GWAS or original GWAS of two compared traits"	
	dat12 = read_ts(snp2gene, trait="DDx")
	dat12 = list(dat12, read_ts(snp2gene_t1, t1), read_ts(snp2gene_t2, t2)) %>% reduce(left_join)		
	dat12 = dat12[order(dat12$FDR_DDx, decreasing=FALSE), ]	
	format_header(wb, "T12.magma_tissue", c(1, 8), h12, header_style)
	format_header(wb, "T12.magma_tissue", c(2, 8), subh12, subheader_style, heights=13*2)
	format_main(wb, "T12.magma_tissue", dat12, title_style, main_style)
	addStyle(wb, "T12.magma_tissue", other_style, rows = which(dat12$FDR_DDx<0.05)+3, cols = 1, gridExpand = TRUE)
	
	# T13.magma_cell_type
	read_cell <- function(fuma_dir, trait)
	{
		dat = fread(paste0(fuma_dir, "/magma_celltype_step1.txt"))
		dat = dat[, c("Dataset", "Cell_type", "P", "P.adj")]
		colnames(dat) = c("Dataset", "Cell_type", paste0(c("P_", "FDR_"), trait))
		return(dat)
	}
	h13 = "Table 13. Significant cell types based on DDx-derived GWAS using MAGMA per dataset cell type analysis"
	subh13 = "Significantly associated cell types at FDR threshold 0.05 were highlighted. FDR for each cell type was estimated by MAGMA based on DDx-derived GWAS or original GWAS of two compared traits; Dataset=scRNA-seq dataset with analyzed cell types."		
	dat13 = read_cell(celltype, trait="DDx")
	dat13 = list(dat13, read_cell(celltype_t1, t1), read_cell(celltype_t2, t2)) %>% reduce(left_join)		
	dat13 = dat13[order(dat13$FDR_DDx, decreasing=FALSE), ]	
	format_header(wb, "T13.magma_cell_type", c(1, ncol(dat13)), h13, header_style)
	format_header(wb, "T13.magma_cell_type", c(2, ncol(dat13)), subh13, subheader_style, heights=13*3)
	format_main(wb, "T13.magma_cell_type", dat13, title_style, main_style)
	addStyle(wb, "T13.magma_cell_type", other_style, rows = which(dat13$FDR_DDx<0.05)+3, cols = 1, gridExpand = TRUE)
	
	# F1.1
	figure1 = paste0("02_ddx/04_FUMA_Results/Image/05_Cell_Type/step1_",pheno,".png")
	if(file.exists(figure1)){
		fh1 = "F1.1 Significant cell types based on DDx-derived GWAS using MAGMA per dataset cell type analysis"	
		format_header(wb, "F1.magma_cell_type", c(1, 14), fh1, header_style)
		insertImage(wb, "F1.magma_cell_type", figure1, startRow = 2, startCol = 1, width = 5, height = 2.5)
	} else{
		removeWorksheet(wb, "F1.magma_cell_type")
	}

	save_path = paste0(save_dir, "/", pheno, ".xlsx")
	saveWorkbook(wb, file = save_path, overwrite = TRUE)
	
	lapply(c(snp2gene, snp2gene_t1, snp2gene_t2, celltype, celltype_t1, celltype_t2), unlink, recursive = TRUE)
	print(paste(pheno, "done!"))
	
	return(save_path)
}

# covidhgi = c("A2", "B2", "C2")
# respiratory = c("Asthma", "ILD", "IPF", "COPD", "Pneumonia.meta", "Pneumonia.b.meta", "Pneumonia.v.meta", "Flu.meta")

# for(t1 in covidhgi)
# {
	# for(t2 in respiratory){ 
		# save_dir = paste0("02_ddx/11_Excel_Summary/Tables/", t1, "/respiratory"); dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
		# DDx_Format_Excel(t1, t2, save_dir)
	# }
# }

