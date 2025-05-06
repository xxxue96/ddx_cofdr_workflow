# Cofdr workflow for COVID-related gwas summary statistics
source("01_data/00_Process_Sumstats.r")
source("03_cofdr/01_02_Run_Cofdr.r")
source("04_script/04_SNP_Clump.r")
source("03_cofdr/02_01_Run_gwas_pw.r")
source("03_cofdr/02_02_gwas_pw_QC.r")
source("03_cofdr/03_01_SNP_Split.r")
source("03_cofdr/03_02_Run_Hyprcoloc.r")
source("03_cofdr/04_SNP_Compare.r")
source("03_cofdr/05_01_Run_Gene_Based_Cofdr_MAGMA.r")
source("03_cofdr/05_02_Run_Gene_Based_Cofdr_PrediXcan.r")
source("03_cofdr/05_03_Run_Gene_Based_Cofdr_PWAS.r")
source("03_cofdr/05_04_Gene_Compare.r")
source("03_cofdr/Cofdr_Workflow_params.r")

# cofdr snp-level analysis
cofdr_wrap_up <- function(t1, t2)
{
	# basic info
	path1 = paste0("01_data/00_Standardise_GWAS/", t1, ".txt.gz")
	path2 = paste0("01_data/00_Standardise_GWAS/", t2, ".txt.gz")
	pheno = paste0(t1, "_", t2)
	traits = c(t1,t2); gwasFilePaths = list(path1,path2)
	merged_path = merge_sumstats(gwasFilePaths, traits)
	overlap_cor = ldsc_intercept_src[t1, t2]
		
	# run decor+cofdr (only applies to overlap_cor<=0.3), extract sig.snp and prepare for FUMA input
	if(overlap_cor<=0.3){
		cor.mat = as.matrix(ldsc_intercept_src[traits, traits])
		cofdr_gwas = Run_Cofdr(traits, gwasFilePaths, merged_path, cor.mat, save_path=paste0("03_cofdr/01_Cofdr_Result/", pheno, ".txt.gz"), sig_save_path=paste0("03_cofdr/02_Cofdr_SigSNP/", pheno, "_", TTthres, "TT.txt.gz"), fuma_save_path=paste0("03_cofdr/05_01_Cofdr_FUMA_Input/", pheno, ".txt.gz"))
		cofdr_snp = cofdr_gwas$sig
	} else{
		print(paste("cofdr not ruuning because of high overlap_cor = ", overlap_cor))
	}
	
	# clumping
	tryCatch({
		clump_gwas = SNP_Clump(path=cofdr_snp, save_path=paste0("03_cofdr/03_Cofdr_Clump/", pheno, ".txt.gz"), clump_kb, clump_r2, clump_p, pop, bfile, plink_bin)
	}, error=function(e){cat("SNP_Clump Warning :",conditionMessage(e), "\n")})
	
	# run gwas-pw
	gwas_pw_res = Run_gwas_pw(t1, t2, path1, path2, save_path=paste0("03_cofdr/06_gwas_pw_Result/", pheno), ldetect_data, overlap_cor)	
	gwas_pw_snp = gwas_pw_QC(gwas_pw_res, chunkPPA3, snpPPA3, save_path=paste0("03_cofdr/07_gwas_pw_SigSNP/", pheno))

	# run hyprcoloc
	## binary.outcomes: continous(0) and binary(1)
	gwas_split_paths = lapply(gwasFilePaths, SNP_Split, region_path=ldetect_data, return_data=TRUE, save_path = tempdir())
	hyprcoloc_res = Run_HyPrColoc(traits, gwas_list=gwas_split_paths, binary.outcomes, save_path=paste0("03_cofdr/08_HyPrColoc_Result/", pheno, ".txt.gz"))

	# compare cofdr, gwas-pw and hyprcoloc 
	tryCatch({
		sigSNP_comparison = SNP_Compare(snp_path=list(cofdr_snp, gwas_pw_snp$sig, hyprcoloc_res), annot_sumstats=merged_path, save_path=paste0("03_cofdr/09_sigSNP_Comparison/", pheno, ".txt.gz"))
	}, error=function(e){cat("SNP_Compare Warning :",conditionMessage(e), "\n")})
}

# cofdr gene-level analysis, when gene-level results for each trait are not available
cofdr_wrap_up_gene <- function(t1, t2)
{
	# basic info
	path1 = paste0("01_data/00_Standardise_GWAS/", t1, ".txt.gz")
	path2 = paste0("01_data/00_Standardise_GWAS/", t2, ".txt.gz")
	pheno = paste0(t1, "_", t2)
	traits = c(t1,t2); gwasFilePaths = list(path1,path2)
	merged_path = merge_sumstats(gwasFilePaths, traits)
	overlap_cor = ldsc_intercept_src[t1, t2]
	
	if(overlap_cor<=0.3){
		cor.mat = as.matrix(ldsc_intercept_src[traits, traits])
		
		twas_dir = paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno); dir.create(twas_dir, showWarnings = FALSE, recursive = TRUE)
		magma_dir=paste0(twas_dir, "/01_MAGMA"); dir.create(magma_dir, showWarnings = FALSE, recursive = TRUE)
		pwas_dir = paste0(twas_dir, "/04_pwas"); dir.create(pwas_dir, showWarnings = FALSE, recursive = TRUE)
		
		# decor + TWAS/spTWAS + cofdr (GRCh38-based mashr model)
		## no imputation due to time limit, so not all model snps available in the decorrelated gwas https://github.com/hakyimlab/MetaXcan/wiki/Best-practices-for-integrating-GWAS-and-GTEX-v8-transcriptome-prediction-models (see Conclusion)	
		cofdr_gene_magma = Run_Gene_Based_Cofdr_MAGMA(traits, gwasFilePaths, merged_path, cor.mat, magma_dir, save_dir=twas_dir)
		cofdr_gene_metaxcan = Run_Gene_Based_Cofdr_MetaXcan(traits, gwasFilePaths, merged_path, cor.mat, save_dir=twas_dir)	
		
		# decor + PWAS + cofdr
		## no imputation due to time limit, so not all LDRef snps available in the decorrelated gwas, and some genes skipped		
		cofdr_gene_pwas = Run_Gene_Based_Cofdr_PWAS(traits, gwasFilePaths, merged_path, cor.mat, pwas_dir, save_dir=twas_dir)
		
		# compare TWAS, spTWAS and PWAS
		tryCatch({
			sigGene_comparison = Gene_Compare(traits, magma_path=cofdr_gene_magma, twas_path=cofdr_gene_metaxcan$eqtl, sptwas_path=cofdr_gene_metaxcan$sqtl, pwas_path=cofdr_gene_pwas, save_path=paste0(twas_dir,"/",pheno,"_geneCompare.txt.gz"))
		}, error=function(e){cat("Gene_Compare Warning :",conditionMessage(e), "\n")})
			
		file.remove(paste0(twas_dir,"/01_harmonized_gwas/",traits, ".txt.gz"))
	} else{
		print(paste("gene-based cofdr not ruuning because of high overlap_cor = ", overlap_cor))
	}
}

# cofdr gene-level analysis, when gene-level results for each trait are available in src="01_data/03_Gene_Result" (applied in a special situation when overlap_cor=0, i.e. cor.mat=NULL)
cofdr_wrap_up_gene2 <- function(t1, t2)
{
	traits = c(t1,t2)
	pheno = paste0(traits, collapse = "_")
	save_dir = twas_dir = paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno)
	magma_dir = paste0(twas_dir, "/01_MAGMA"); dir.create(magma_dir, showWarnings = FALSE, recursive = TRUE)
	pwas_dir = paste0(twas_dir, "/04_pwas"); dir.create(pwas_dir, showWarnings = FALSE, recursive = TRUE)
	src = "01_data/03_Gene_Result"
	
	magmaPaths = as.list(paste0(src, "/", traits, "/01_MAGMA","/",traits,".0UP.0DOWN/Entrze/",traits,".0UP.0DOWN.genes.out"))
	cofdr_gene_magma = run_magma_cofdr(traits, magmaPaths, magma_dir, save_dir)

	predixcan_eqtl = c(list.files(paste0(src, "/", t1, "/02_spredixcan/eqtl/",t1), pattern=".csv", full.names=TRUE), list.files(paste0(src, "/", t2,"/02_spredixcan/eqtl/",t2), pattern=".csv", full.names=TRUE))
	predixcan_sqtl = c(list.files(paste0(src, "/", t1,"/02_spredixcan/sqtl/",t1), pattern=".csv", full.names=TRUE), list.files(paste0(src, "/", t2,"/02_spredixcan/sqtl/",t2), pattern=".csv", full.names=TRUE))
	multixcan_eqtl = c(list.files(path=paste0(src, "/", t1, "/03_smultixcan/eqtl"), pattern=".txt$", full.names=TRUE), list.files(path=paste0(src, "/", t2, "/03_smultixcan/eqtl"), pattern=".txt$", full.names=TRUE))
	multixcan_sqtl = c(list.files(path=paste0(src, "/", t1, "/03_smultixcan/sqtl"), pattern=".txt$", full.names=TRUE), list.files(path=paste0(src, "/", t2, "/03_smultixcan/sqtl"), pattern=".txt$", full.names=TRUE))
	metaxcanPaths = list(predixcan_eqtl=predixcan_eqtl, predixcan_sqtl=predixcan_sqtl, multixcan_eqtl=multixcan_eqtl, multixcan_sqtl=multixcan_sqtl)
	dir.create(paste0(twas_dir, "/02_spredixcan/eqtl"), showWarnings = FALSE, recursive = TRUE)
	dir.create(paste0(twas_dir, "/02_spredixcan/sqtl"), showWarnings = FALSE, recursive = TRUE)
	dir.create(paste0(twas_dir, "/03_smultixcan/eqtl"), showWarnings = FALSE, recursive = TRUE)
	dir.create(paste0(twas_dir, "/03_smultixcan/sqtl"), showWarnings = FALSE, recursive = TRUE)	
	cofdr_gene_metaxcan = run_metaxcan_cofdr(traits, metaxcanPaths, save_dir)

	pwas_stat <- function(i)
	{
		pheno=traits[i]
		res = paste0(src, "/", pheno, "/04_pwas/", pheno, "_chr", 1:22, ".out")
		res = do.call(rbind.data.frame, lapply(res, fread))
		res = res[, c("ID", "PWAS.Z")]		
		colnames(res) = c("SNP", paste0("Z_",pheno))
			
		return(res)
	}	
	pwasRes = lapply(1:length(traits), pwas_stat)
	cofdr_gene_pwas = run_pwas_cofdr(traits, pwasRes, pwas_dir, save_dir)
	
	tryCatch({
		sigGene_comparison = Gene_Compare(traits, magma_path=cofdr_gene_magma, twas_path=cofdr_gene_metaxcan$eqtl, sptwas_path=cofdr_gene_metaxcan$sqtl, pwas_path=cofdr_gene_pwas, save_path=paste0(twas_dir,"/",pheno,"_geneCompare.txt.gz"))
	}, error=function(e){cat("Gene_Compare Warning :",conditionMessage(e), "\n")})
}

# args <- commandArgs(TRUE)
# t1 = args[1]
# t2 = args[2]

# if(!is.na(t1) & !is.na(t2)){
	# cofdr_wrap_up(t1, t2)
	# # cofdr_wrap_up_gene(t1, t2)
	# cofdr_wrap_up_gene2(t1, t2)
# }

