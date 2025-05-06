# colocalized pGenes (i.e. map snps to genes based on pqtl)
library(data.table)
library(dplyr)
setDTthreads(40)
source("01_data/00_Process_Sumstats.r")
source("03_cofdr/01_00_Zmat_Decor.r")
source("03_cofdr/01_02_Run_Cofdr.r")
source("03_cofdr/Cofdr_Workflow_params.r")

# step1: decorrelate gwas summary statistics 
# @decor_path: directory to save decorrelated gwas (A1 as effect allele)
decor_snp4pwas <- function(traits, gwasFilePaths, merged_path=NULL, cor.mat=NULL, decor_path=tempdir())
{
	common_snp = merge_sumstats(gwasFilePaths, traits, merged_path)[, c("SNP", "A2", "A1", paste0("Z_", traits)), with=FALSE]
	colnames(common_snp)[2:3] = c("A1", "A2")
	zmat = as.matrix(common_snp[, c("SNP", grep("^Z_", colnames(common_snp), value=TRUE)), with=FALSE], rownames="SNP")
	if(!is.null(cor.mat)) {zmat = t(decor(t(zmat), cor.mat))}

	## save decorrelated gwas with SNP, A1(effect_allele), A2, Z
	common_snp = cbind(as.data.frame(common_snp[, c("SNP", "A1", "A2")]), as.data.frame(zmat))		
	create_decor_gwas_paths <- function(pheno)
	{
		decor_gwas = common_snp[, c("SNP", "A1", "A2", paste0("Z_", pheno))]
		colnames(decor_gwas) = c("SNP", "A1", "A2", "Z")
				
		outPath = paste0(decor_path, "/", pheno, ".txt")
		fwrite(decor_gwas, file=outPath, sep=" ", quote=FALSE)
		return(outPath)
	}
	decor_gwas_paths = lapply(traits, create_decor_gwas_paths)
	
	return(decor_gwas_paths)
}


# step2: run pwas 
# @decor_gwas_paths: decorrelated gwas paths
# @pwas_dir: directory to save pwas results
create_pwas_gene_results <- function(traits, decor_gwas_paths, pwas_dir)
{
	run_pwas <- function(i)
	{
		pheno=traits[i]; INPUT=decor_gwas_paths[[i]]; OUTPUT=pwas_dir		
		cmd = paste("bash 04_script/07_Run_PWAS.bash", pheno, INPUT, OUTPUT)
		system(cmd)
		
		res = paste0(pwas_dir, "/", pheno, "_chr", 1:22, ".out")
		res = do.call(rbind.data.frame, lapply(res, fread))
		res = res[, c("ID", "PWAS.Z")]
		colnames(res) = c("SNP",  paste0("Z_",pheno))
		
		return(res)
	}
	pwasRes = lapply(1:length(traits), run_pwas)
	
	return(pwasRes)
}


# step3: run gene-based cofdr with zmat obtained from step2
# @pwas_dir: directory to save PWAS results
# @save_dir: directory to save significant cofdr pGenes
run_pwas_cofdr <- function(traits, pwasRes, pwas_dir, save_dir)
{	
	pheno = paste(traits, collapse="_")
	
	## run cofdr for PWAS gene zscores, remove NA/Inf for zscore cols
	merged_path = pwasRes %>% reduce(inner_join, by = "SNP")
	merged_path = merged_path[is.finite(rowSums(merged_path[, -"SNP"])), ]
	tryCatch({cofdr_pqtl = Run_Cofdr(traits, gwasFilePaths=NULL, merged_path, cor.mat=NULL, save_path=paste0(pwas_dir, "/", pheno, ".txt.gz"), keep_Zstat=TRUE)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
		
	if(exists("cofdr_pqtl")){
		## filter and annotate sig.cofdr results
		cofdr_pqtl = cofdr_pqtl$allsnp
		cofdr_pqtl = fread(cofdr_pqtl)[, c("SNP", "TT")]
		sig_cofdr_pqtl = cofdr_pqtl[cofdr_pqtl$TT>=TTthres_pqtl, ]
		colnames(sig_cofdr_pqtl)[1] <- "hgnc_symbol"
		
		## only keep the gene with largest TT
		sig_cofdr_pqtl = sig_cofdr_pqtl %>% group_by(hgnc_symbol) %>% filter(TT == max(TT, na.rm=TRUE))	%>% as.data.table	
	} else{
		## not enough genes, and thus failed in running cofdr
		header = c("SNP", "TT")
		sig_cofdr_pqtl = setNames(data.frame(matrix(ncol = length(header), nrow = 0)), header)
	}
		
	sig_save_path = paste0(save_dir,"/",pheno,"_sig.pqtl.txt.gz")
	fwrite(sig_cofdr_pqtl, file=sig_save_path, sep="\t")
	return(sig_save_path)
}


# wrap-up functions
Run_Gene_Based_Cofdr_PWAS <- function(traits, gwasFilePaths, merged_path, cor.mat=NULL, pwas_dir, save_dir)
{	
	decor_gwas_paths = decor_snp4pwas(traits, gwasFilePaths, merged_path, cor.mat, decor_path=tempdir())
	pwasRes = create_pwas_gene_results(traits, decor_gwas_paths, pwas_dir)	
	gene_based_cofdr = run_pwas_cofdr(traits, pwasRes, pwas_dir, save_dir)
	
	return(gene_based_cofdr)
}