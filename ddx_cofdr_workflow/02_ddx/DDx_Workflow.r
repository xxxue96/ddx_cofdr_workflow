source("01_data/00_Process_Sumstats.r")
source("02_ddx/02_Run_DDx.r")
source("04_script/04_SNP_Clump.r")
source("02_ddx/03_02_Run_CCGWAS.r")
source("02_ddx/05_01_Run_DDx_Gene_MAGMA.r")
source("02_ddx/05_02_Run_DDx_Gene_PrediXcan.r")
source("02_ddx/05_03_Run_DDx_Gene_PWAS.r")
source("02_ddx/05_04_Gene_Compare.r")
source("02_ddx/DDx_Workflow_params.r") ## default parameters can be changed

# ddx snp-level analysis
ddx_wrap_up <- function(t1, t2)
{
	# basic info
	path1 = paste0("01_data/00_Standardise_GWAS/", t1, ".txt.gz")
	path2 = paste0("01_data/00_Standardise_GWAS/", t2, ".txt.gz")
	pheno = paste0(t1, "_", t2)
	merged_path = merge_sumstats(gwasFilePaths=list(path1,path2), traits=c(t1,t2))
	ldsc_intercept = as.numeric(ldsc_intercept_src[t1, t2])
		
	# run ddx and extract sig.snp
	ddx_snp = paste0("02_ddx/02_DDx_SigSNP/", pheno, "_", pthres, ".txt.gz")
	ddx_gwas = Run_DDx(t1, t2, path1, path2, merged_path, ldsc_intercept, save_path=paste0("02_ddx/01_DDx_Results/", pheno, ".txt.gz"), sig_save_path=ddx_snp)

	# clumping
	tryCatch({
		clump_gwas = SNP_Clump(path=ddx_snp, save_path=paste0("02_ddx/03_DDx_Clump/", pheno, ".txt.gz"), clump_kb, clump_r2, clump_p, pop, bfile, plink_bin)
	}, error=function(e){cat("SNP_Clump Warning:",conditionMessage(e), "\n")})
	
	# run ccgwas (subtype_data=TRUE when comparing subtypes of a disorder)
	ccgwas_snp = paste0("02_ddx/06_CCGWAS_SigSNP/", pheno, ".txt.gz")
	ccgwas_res = Run_CCGWAS(t1, t2, path1, path2, save_path=paste0("02_ddx/05_CCGWAS_Results/", pheno), subtype_data=FALSE, sig_save_path=ccgwas_snp)
	
	# compare ddx and ccgwas significant snp	
	tryCatch({
		ddx_vs_ccgwas = merge_sumstats(gwasFilePaths=list(ddx_snp, ccgwas_snp), traits=NA, merged_path=NULL, save_path=paste0("02_ddx/07_SNP_Compare/", pheno, ".txt.gz"))
	}, error=function(e){cat("ddx_vs_ccgwas Warning:",conditionMessage(e), "\n")})
	
	tryCatch({
		ddx_vs_ccgwas_clump = SNP_Clump(path=ddx_vs_ccgwas, save_path=paste0("02_ddx/07_SNP_Compare/", pheno, ".clump.txt.gz"), clump_kb, clump_r2, clump_p, pop, bfile, plink_bin)
	}, error=function(e){cat("ddx_vs_ccgwas_clump Warning:",conditionMessage(e), "\n")})
}


# ddx gene-level analysis
# @ddx_name: set a pseudo-name for the DDx-derived trait sumstats 
ddx_wrap_up_gene <- function(t1, t2, ddx_name=NULL)
{	
	# basic info	
	traits = c(t1, t2)
	if(is.null(ddx_name)){ pheno = paste0(traits, collapse="_") } else{ pheno = ddx_name }
	
	# ddx sumstats reformat
	ddx_gwas = paste0("02_ddx/01_DDx_Results/", pheno, ".txt.gz")
	twas_dir_ddx = paste0("02_ddx/08_DDx_Gene_Result/", pheno)
	magma_dir_ddx = paste0(twas_dir_ddx, "/01_MAGMA")
	pwas_dir_ddx = paste0(twas_dir_ddx, "/04_pwas")
	sapply(list(twas_dir_ddx, magma_dir_ddx, pwas_dir_ddx), dir.create, showWarnings = FALSE)
	
	# magma
	ddx_gene_magma = Run_DDx_Gene_MAGMA(decor_gwas_paths=ddx_gwas, magma_dir=magma_dir_ddx, save_dir=twas_dir_ddx)
	
	# predixcan
	ddx_gene_metaxcan = Run_DDx_Gene_PrediXcan(traits=pheno, decor_gwas_paths=ddx_gwas, metaxcan_dir=twas_dir_ddx, save_dir=twas_dir_ddx)
	
	# pwas
	ddx_gene_pwas = Run_DDx_Gene_PWAS(traits=pheno, decor_gwas_paths=ddx_gwas, pwas_dir=pwas_dir_ddx, save_dir=twas_dir_ddx)
	
	# compare ddx_gene_magma/metaxcan/pwas
	magma_gene = fread(ddx_gene_magma) %>% rename_at(vars(-hgnc_symbol), ~ paste0(., '_DDx')) %>% mutate_at("hgnc_symbol", as.character)
	twas_gene = fread(ddx_gene_metaxcan$eqtl)%>% rename(P=pvalue, Z=zscore) %>% rename_at(vars(-hgnc_symbol, -tissue), ~ paste0(., '_DDx'))%>% mutate_at(c("hgnc_symbol","tissue"), as.character) %>% relocate(hgnc_symbol, tissue)	
	sptwas_gene = fread(ddx_gene_metaxcan$sqtl) %>% rename(P=pvalue, Z=zscore) %>% dplyr::select(-ensemble_id) %>% rename_at(vars(-hgnc_symbol, -tissue, -sqtl_hg38, -sqtl_hg19), ~ paste0(., '_DDx'))%>% mutate_at(c("hgnc_symbol","tissue","sqtl_hg38","sqtl_hg19"), as.character) %>% relocate(hgnc_symbol, tissue)	
	pwas_gene = fread(ddx_gene_pwas) %>% rename(P=PWAS.P, Z=PWAS.Z) %>% dplyr::select(-FILE) %>% rename_at(vars(-hgnc_symbol), ~ paste0(., '_DDx')) %>% mutate_at("hgnc_symbol", as.character)
	sigGene_comparison_ddx = Gene_Compare(traits, magma_gene, twas_gene, sptwas_gene, pwas_gene, save_path=paste0("02_ddx/08_DDx_Gene_Result/",pheno,"/",pheno,"_geneCompare.txt.gz"))

	# when CC-GWAS results available
	if(file.exists(paste0("02_ddx/05_CCGWAS_Results/", pheno, ".results.gz")))
	{
		# ccgwas sumstats reformat
		reformat_ccgwas_res <- function(t1, t2, ccgwas_res_path)
		{
			dat = fread(ccgwas_res_path)[, c("SNP", "CHR", "BP", "NEA", "EA", "OLS_beta", "OLS_se", "OLS_pval")]
			colnames(dat) = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P")
			dat$Z = dat$BETA / dat$SE
			dat$N = as.numeric(sample_size_src[t1, "N"] + sample_size_src[t2, "N"])
					
			twas_dir_ccgwas = paste0("02_ddx/09_CCGWAS_Gene_Result/", t1, "_", t2, "/OLS")
			save_dir = paste0(twas_dir_ccgwas, "/00_tmp")
			magma_dir_ccgwas = paste0(twas_dir_ccgwas, "/01_MAGMA")
			pwas_dir_ccgwas = paste0(twas_dir_ccgwas, "/04_pwas")
			sapply(list(twas_dir_ccgwas, magma_dir_ccgwas, pwas_dir_ccgwas, save_dir), dir.create, showWarnings = FALSE, recursive = TRUE)
			
			save_path = paste0(save_dir, "/", t1, "_", t2, ".txt.gz")
			fwrite(dat, file=save_path, sep="\t", quote=FALSE, na=NA)
		
			return(list(save_path=save_path, twas_dir=twas_dir_ccgwas, magma_dir=magma_dir_ccgwas, pwas_dir=pwas_dir_ccgwas))
		}
		ccgwas_ols = reformat_ccgwas_res(t1, t2, ccgwas_res_path=paste0("02_ddx/05_CCGWAS_Results/", pheno, ".results.gz"))
		
		# magma
		ccgwas_ols_gene_magma = Run_DDx_Gene_MAGMA(decor_gwas_paths=ccgwas_ols$save_path, magma_dir=ccgwas_ols$magma_dir, save_dir=ccgwas_ols$twas_dir)
		
		# predixcan
		ccgwas_ols_gene_metaxcan = Run_DDx_Gene_PrediXcan(traits=pheno, decor_gwas_paths=ccgwas_ols$save_path, metaxcan_dir=ccgwas_ols$twas_dir, save_dir=ccgwas_ols$twas_dir)
		
		# pwas	
		ccgwas_ols_gene_pwas = Run_DDx_Gene_PWAS(traits=pheno, decor_gwas_paths=ccgwas_ols$save_path, pwas_dir=ccgwas_ols$pwas_dir, save_dir=ccgwas_ols$twas_dir)
			
		# compare ccgwas_ols_gene_magma/metaxcan/pwas
		magma_gene = fread(ccgwas_ols_gene_magma) %>% rename_at(vars(-hgnc_symbol), ~ paste0(., '_DDx')) %>% mutate_at("hgnc_symbol", as.character)
		twas_gene = fread(ccgwas_ols_gene_metaxcan$eqtl)%>% rename(P=pvalue, Z=zscore) %>% rename_at(vars(-hgnc_symbol, -tissue), ~ paste0(., '_DDx'))%>% mutate_at(c("hgnc_symbol","tissue"), as.character) %>% relocate(hgnc_symbol, tissue)	
		sptwas_gene = fread(ccgwas_ols_gene_metaxcan$sqtl) %>% rename(P=pvalue, Z=zscore) %>% dplyr::select(-ensemble_id) %>% rename_at(vars(-hgnc_symbol, -tissue, -sqtl_hg38, -sqtl_hg19), ~ paste0(., '_DDx'))%>% mutate_at(c("hgnc_symbol","tissue","sqtl_hg38","sqtl_hg19"), as.character) %>% relocate(hgnc_symbol, tissue)	
		pwas_gene = fread(ccgwas_ols_gene_pwas) %>% rename(P=PWAS.P, Z=PWAS.Z) %>% dplyr::select(-FILE) %>% rename_at(vars(-hgnc_symbol), ~ paste0(., '_DDx')) %>% mutate_at("hgnc_symbol", as.character)
		sigGene_comparison_ccgwas_ols = Gene_Compare(traits, magma_gene, twas_gene, sptwas_gene, pwas_gene, save_path=paste0("02_ddx/09_CCGWAS_Gene_Result/",pheno,"/",pheno,"_geneCompare.txt.gz"))
		tmp = fread(sigGene_comparison_ccgwas_ols)
		colnames(tmp) = gsub("_DDx", "_OLS", colnames(tmp))
		fwrite(tmp, file=paste0("02_ddx/09_CCGWAS_Gene_Result/",pheno,"/OLS/",pheno,"_geneCompare.txt.gz"), sep="\t")
		
		# first compare gene results based on DDx and CCGWAS_OLS effect sizes; then compare magma_gene, twas_gene, sptwas_gene, and pwas_gene
		magma_gene = merge_magma(ddx_gene_magma, ccgwas_ols_gene_magma)
		twas_gene = merge_twas(ddx_gene_metaxcan$eqtl, ccgwas_ols_gene_metaxcan$eqtl)
		sptwas_gene = merge_sptwas(ddx_gene_metaxcan$sqtl, ccgwas_ols_gene_metaxcan$sqtl)
		pwas_gene = merge_pwas(ddx_gene_pwas, ccgwas_ols_gene_pwas)	
		gene_validate = Gene_Compare(traits, magma_gene, twas_gene, sptwas_gene, pwas_gene, save_path=paste0("02_ddx/10_DDx_CCGWAS_Gene_Compare/",pheno,"_geneCompare.txt.gz"))
	}
	
	# functional annotation of identified genes
	## FUMA：GSEA; magma gene-set, tissue and cell type enrichment
	
	# remove temporary files
	tryCatch({
		file.remove(c(ccgwas_ols$save_path, ccgwas_exact$save_path, paste0(c(twas_dir_ddx, ccgwas_ols$twas_dir, ccgwas_exact$twas_dir),"/01_harmonized_gwas/",pheno,".txt.gz")))
	}, error=function(e){cat("FILE REMOVE Warning:",conditionMessage(e), "\n")})
	
	return(sigGene_comparison_ddx)
}


# args <- commandArgs(TRUE)
# t1 = args[1]
# t2 = args[2]

# if(!is.na(t1) & !is.na(t2)){
	# ddx_wrap_up(t1, t2)
	# ddx_wrap_up_gene(t1, t2)
# }