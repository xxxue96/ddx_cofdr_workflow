library(openxlsx)
library(data.table)
library(dplyr)
source("Excel_Format_Params.r")
source("02_ddx/mtCOJO_Workflow_params.r")

# @mtcojo_name: set a pseudo-name for the mtCOJO-derived trait sumstats 
mtcojo_Format_Excel <- function(t1, t2phenos, mtcojo_name=NULL, save_dir)
{
	traits = c(t1, t2phenos)
	if(is.null(mtcojo_name)){ pheno = paste0(traits, collapse="_") } else{ pheno = mtcojo_name }
	print(paste(pheno, "start!"))
	snp2gene = paste0("02_ddx/04_FUMA_Results/SNP2GENE/",pheno)
	snp2gene_ann = paste0("01_data/04_FUMA_Results/SNP2GENE/",traits)
	celltype = paste0("02_ddx/04_FUMA_Results/Cell_Type/",pheno)
	celltype_ann = paste0("01_data/04_FUMA_Results/Cell_Type/",traits)
	image_dir = "04_FUMA_Results/Image/05_Cell_Type"	
	lapply(c(snp2gene, snp2gene_ann, celltype, celltype_ann), function(file_path){unzip(paste0(file_path,".zip"), exdir=file_path)})
	
	# create workbook
	wb <- createWorkbook()
	addWorksheet(wb, "T1.sig_snps")
	addWorksheet(wb, "T2.risk_loci")
	addWorksheet(wb, "T3.GWAS_catalog")
	addWorksheet(wb, "T4.candidate_snps")
	addWorksheet(wb, "T5.annov_enrichment")
	addWorksheet(wb, "T6.mapping_genes")
	addWorksheet(wb, "T7.gene_mtcojo")	
	addWorksheet(wb, "T8.magma_gene_set")
	addWorksheet(wb, "T9.magma_tissue")
	addWorksheet(wb, "T10.magma_cell_type")
	addWorksheet(wb, "F1.magma_cell_type")

	# T1.sig_snps
	h1 = paste("Table 1. Significant SNPs detected by mtCOJO at p-value threshold", pthres)
	subh1 = "A1: non-effect allele; A2: effect allele"
	dat1 = fread(paste0("02_ddx/02_DDx_SigSNP/", pheno, "_", pthres, ".txt.gz"))
	dat1 = dat1[, colnames(dat1)[!grepl("N_", colnames(dat1))], with=FALSE]
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
	
	# T7.gene_mtcojo
	h7 = "Table 7. Results summary of MAGMA/TWAS/spTWAS/PWAS-mtCOJO"
	subh7 = "Significant genes associated with conditional traits were identified by integrating gene genomic positions, mRNA expression, splicing and plasma protein data with mtCOJO-derived GWAS using MAGMA gene analysis, PrediXcan or FUSION. Overlapping results of at least three approaches of MAGMA/TWAS/spTWAS/PWAS-mtCOJO were highlighted (total_detect>=3). P_magma=the gene p-value estimated by MAGMA based on mtCOJO-derived GWAS; FDR_magma=adjusted P_magma; P_twas=the gene p-value estimated by TWAS based on mtCOJO-derived GWAS; tissue_twas=the eQTL model of every single tissue used to obtain gene-based z-scores in TWAS analysis based on mtCOJO-derived GWAS; FDR_twas=adjusted P_twas; zscore_twas=gene-based z-scores in TWAS analysis; P_sptwas=the gene p-value estimated by spTWAS based on mtCOJO-derived GWAS; tissue_sptwas=the sQTL model of every single tissue used to obtain gene-based z-scores in spTWAS analysis based on mtCOJO-derived GWAS; FDR_sptwas=adjusted P_sptwas; zscore_sptwas=gene-based z-scores in spTWAS analysis; P_pwas=the gene p-value estimated by PWAS based on mtCOJO-derived GWAS; FDR_pwas=adjusted P_pwas; zscore_pwas=gene-based z-scores in PWAS analysis"
	dat7 = fread(paste0("02_ddx/08_DDx_Gene_Result/", pheno, "/", pheno, "_geneCompare.txt.gz"))
	dat7 = dat7[order(dat7$total_detect, decreasing=TRUE), ]
	colnames(dat7)[1] = "Gene"
	format_header(wb, "T7.gene_mtcojo", c(1, ncol(dat7)), h7, header_style)
	format_header(wb, "T7.gene_mtcojo", c(2, ncol(dat7)), subh7, subheader_style, heights=13*1)
	format_main(wb, "T7.gene_mtcojo", dat7, title_style, main_style)
	addStyle(wb, "T7.gene_mtcojo", other_style, rows = which(dat7$total_detect>=3)+3, cols = 1, gridExpand = TRUE)
		
	# T8.magma_gene_set
	read_gene_set <- function(fuma_dir, trait)
	{
		dat = tryCatch(
		{
			dat = fread(paste0(fuma_dir, "/magma.gsa.out"), skip="# CONDITIONED_INTERNAL = gene size, gene density, sample size, inverse mac, log(gene size), log(gene density), log(sample size), log(inverse mac)")
		}, error=function(e){
			return(fread(paste0(fuma_dir, "/magma.gsa.out"), skip="# CONDITIONED_INTERNAL = gene size, gene density, inverse mac, log(gene size), log(gene density), log(inverse mac)"))
		})
				
		dat$FDR = p.adjust(dat$P, method = "fdr")
		dat = dat[, c("FULL_NAME", "P", "FDR")]
		colnames(dat) = c("Gene_set", paste0(c("P_", "FDR_"), trait))
		return(dat)
	}	
	h8 = "Table 8. Gene set analysis based on mtCOJO-derived GWAS using MAGMA"
	subh8 = "Significantly associated gene sets at FDR threshold 0.05 were highlighted. FDR for each gene set was estimated by MAGMA based on mtCOJO-derived GWAS or original GWAS of compared traits"
	dat8 = read_gene_set(snp2gene, trait="mtCOJO")
	dat8 = c(list(dat8), lapply(1:length(traits), function(i){read_gene_set(snp2gene_ann[i], traits[i])})) %>% reduce(left_join)	
	dat8 = dat8[order(dat8$FDR_mtCOJO, decreasing=FALSE), ]
	format_header(wb, "T8.magma_gene_set", c(1, ncol(dat8)), h8, header_style)
	format_header(wb, "T8.magma_gene_set", c(2, ncol(dat8)), subh8, subheader_style, heights=13*1)
	format_main(wb, "T8.magma_gene_set", dat8, title_style, main_style)
	addStyle(wb, "T8.magma_gene_set", other_style, rows = which(dat8$FDR_mtCOJO<0.05)+3, cols = 1, gridExpand = TRUE)

	# T9.magma_tissue
	read_ts <- function(fuma_dir, trait)
	{
		dat = fread(paste0(fuma_dir, "/magma_exp_gtex_v8_ts_avg_log2TPM.gsa.out"), skip="# CONDITIONED_HIDDEN = Average (covar)")
		dat$FDR = p.adjust(dat$P, method = "fdr")
		dat = dat[, c("FULL_NAME", "P", "FDR")]
		colnames(dat) = c("Tissue", paste0(c("P_", "FDR_"), trait))
		return(dat)
	}
	h9 = "Table 9. Tissue enrichment analysis based on mtCOJO-derived GWAS using MAGMA"
	subh9 = "Significantly associated tissues at FDR threshold 0.05 were highlighted. FDR for each tissue was estimated by MAGMA based on mtCOJO-derived GWAS or original GWAS of compared traits"	
	dat9 = read_ts(snp2gene, trait="mtCOJO")
	dat9 = c(list(dat9), lapply(1:length(traits), function(i){read_ts(snp2gene_ann[i], traits[i])})) %>% reduce(left_join)	
	dat9 = dat9[order(dat9$FDR_mtCOJO, decreasing=FALSE), ]	
	format_header(wb, "T9.magma_tissue", c(1, ncol(dat9)), h9, header_style)
	format_header(wb, "T9.magma_tissue", c(2, ncol(dat9)), subh9, subheader_style, heights=13*1)
	format_main(wb, "T9.magma_tissue", dat9, title_style, main_style)
	addStyle(wb, "T9.magma_tissue", other_style, rows = which(dat9$FDR_mtCOJO<0.05)+3, cols = 1, gridExpand = TRUE)

	# T10.magma_cell_type
	read_cell <- function(fuma_dir, trait)
	{
		dat = fread(paste0(fuma_dir, "/magma_celltype_step1.txt"))
		dat = dat[, c("Dataset", "Cell_type", "P", "P.adj")]
		colnames(dat) = c("Dataset", "Cell_type", paste0(c("P_", "FDR_"), trait))
		return(dat)
	}
	h10 = "Table 10. Significant cell types based on mtCOJO-derived GWAS using MAGMA per dataset cell type analysis"
	subh10 = "Significantly associated cell types at FDR threshold 0.05 were highlighted. FDR for each cell type was estimated by MAGMA based on mtCOJO-derived GWAS or original GWAS of compared traits; Dataset=scRNA-seq dataset with analyzed cell types."		
	dat10 = read_cell(celltype, trait="mtCOJO")
	dat10 = c(list(dat10), lapply(1:length(traits), function(i){read_cell(celltype_ann[i], traits[i])})) %>% reduce(left_join)		
	dat10 = dat10[order(dat10$FDR_mtCOJO, decreasing=FALSE), ]	
	format_header(wb, "T10.magma_cell_type", c(1, ncol(dat10)), h10, header_style)
	format_header(wb, "T10.magma_cell_type", c(2, ncol(dat10)), subh10, subheader_style, heights=13*1)
	format_main(wb, "T10.magma_cell_type", dat10, title_style, main_style)
	addStyle(wb, "T10.magma_cell_type", other_style, rows = which(dat10$FDR_mtCOJO<0.05)+3, cols = 1, gridExpand = TRUE)
	
	# F1.1
	figure1 = paste0("02_ddx/04_FUMA_Results/Image/05_Cell_Type/step1_",pheno,".png")
	if(file.exists(figure1)){
		fh1 = "F1.1 Significant cell types for mtCOJO-derived association statistics using MAGMA per dataset cell type analysis"	
		format_header(wb, "F1.magma_cell_type", c(1, 14), fh1, header_style)
		insertImage(wb, "F1.magma_cell_type", figure1, startRow = 2, startCol = 1, width = 5, height = 2.5)
	} else{
		removeWorksheet(wb, "F1.magma_cell_type")
	}
	
	save_path = paste0(save_dir, "/", pheno, ".xlsx")
	saveWorkbook(wb, file = save_path, overwrite = TRUE)

	lapply(c(snp2gene, snp2gene_ann, celltype, celltype_ann), unlink, recursive = TRUE)
	print(paste(pheno, "done!"))
}


# covidhgi = c("A2", "B2", "C2")
# t2phenos = c("Asthma", "ILD", "IPF", "COPD", "Pneumonia.meta", "Pneumonia.b.meta", "Pneumonia.v.meta", "Flu.meta")

# for(t1 in covidhgi)
# {
	# mtcojo_Format_Excel(t1, t2phenos, mtcojo_name=paste0(t1,"_respiratory"), save_dir="13_mtCOJO_Excel_Summary/Tables")
# }
