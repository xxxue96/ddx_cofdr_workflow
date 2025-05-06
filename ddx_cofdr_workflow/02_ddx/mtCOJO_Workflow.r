source("02_ddx/DDx_Workflow.r")
source("02_ddx/04_Run_mtcojo.r")
source("02_ddx/mtCOJO_Workflow_params.r")

# mtcojo snp-level analysis
# @mtcojo_name: set a pseudo-name for the mtCOJO-derived trait sumstats 
mtCOJO_wrap_up <- function(t1, t2phenos, mtcojo_name=NULL)
{
	# basic info
	traits = c(t1, t2phenos)
	if(is.null(mtcojo_name)){ pheno = paste0(traits, collapse="_") } else{ pheno = mtcojo_name }
	gwasFilePaths = paste0("01_data/00_Standardise_GWAS/", traits, ".txt.gz")
	merged_path = merge_sumstats(gwasFilePaths, traits)
	merged_path = merged_path[, c("SNP", "CHR", "BP", "A1", "A2", as.character(outer(X=c("BETA_", "SE_", "P_"), Y=traits, FUN=paste0))), with=FALSE]
	
	# run mtcojo and extract sig.snp
	ddx_snp = paste0("02_ddx/02_DDx_SigSNP/", pheno, "_", pthres, ".txt.gz")
	ddx_gwas = Run_mtcojo(t1, t2phenos, gwasFilePaths, merged_path, save_path=paste0("02_ddx/12_mtCOJO_Results/", pheno), sig_save_path=ddx_snp)
	
	# clumping	
	tryCatch({
		clump_gwas = SNP_Clump(path=ddx_snp, save_path=paste0("02_ddx/03_DDx_Clump/", pheno, ".txt.gz"), clump_kb, clump_r2, clump_p, pop, bfile, plink_bin)
	}, error=function(e){cat("SNP_Clump Warning",conditionMessage(e), "\n")})	
}

# mtcojo gene-level analysis
mtCOJO_wrap_up_gene <- function(t1, t2phenos, mtcojo_name=NULL)
{
	# basic info
	traits = c(t1, t2phenos)
	if(is.null(mtcojo_name)){ pheno = paste0(traits, collapse="_") } else{ pheno = mtcojo_name }
	file.copy(from=paste0("02_ddx/12_mtCOJO_Results/", pheno,".txt.gz"), to=paste0("02_ddx/01_DDx_Results/", pheno,".txt.gz"))
	
	# compare magma/metaxcan/pwas results based on mtcojo-derived sumstats
	sigGene_comparison_ddx = ddx_wrap_up_gene(t1, t2=t2phenos, ddx_name=pheno)
	tmp = fread(sigGene_comparison_ddx) %>% rename_with(~ gsub("DDx", "mtCOJO", .x))
	fwrite(tmp, file=sigGene_comparison_ddx, sep="\t")
}

# covidhgi = c("A2", "B2", "C2")
# t2phenos = c("Asthma", "ILD", "IPF", "COPD", "Pneumonia.meta", "Pneumonia.b.meta", "Pneumonia.v.meta", "Flu.meta")

# for(t1 in covidhgi){
	# mtCOJO_wrap_up(t1, t2phenos, mtcojo_name=paste0(t1,"_respiratory"))
	# mtCOJO_wrap_up_gene(t1, t2phenos, mtcojo_name=paste0(t1,"_respiratory"))
# }