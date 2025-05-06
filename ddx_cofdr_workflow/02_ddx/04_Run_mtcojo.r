# run mtcojo to identify disorder-specific variants https://yanglab.westlake.edu.cn/software/gcta/#mtCOJO
library(data.table)
setDTthreads(40)
source("01_data/00_Process_Sumstats.r")
source("01_data/01_Standardise_GWAS.r")
source("04_script/04_SNP_Clump.r")
source("02_ddx/mtCOJO_Workflow_params.r")

# prepare data input for "--mtcojo-file"
# the sample prevalence and the population prevalence are not specified, the estimate of the SNP-based h2 will be on the observed scale.
cojo_data_prepare <- function(path, save_path=tempdir())
{
	# reformat
	gwas = read_sumstats(path)[, c("SNP", "A2", "A1", "FRQ", "BETA", "SE", "P", "N")]
	colnames(gwas) = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
	
	# save
	trait = gsub(".txt.gz", "", basename(path))	
	save_path = paste0(save_path, "/", trait, ".txt.gz")
	fwrite(gwas, file=save_path, sep="\t", quote=FALSE)
	
	return(c(trait, save_path))	
}

# @gwasFilePaths: The first path is for the target trait (i.e. y), and the remaining paths are for covariate traits.
Run_mtcojo <- function(t1, t2phenos, gwasFilePaths, merged_path=NULL, save_path, sig_save_path)
{
	# basic info
	traits = c(t1, t2phenos)
	merged_path = merge_sumstats(gwasFilePaths, traits, merged_path)
	
	# bfile
	ref_geno = bfile
	
	# mtcojo-file
	gwas = do.call(rbind, lapply(gwasFilePaths, cojo_data_prepare))
	gwas_path = paste0(tempdir(), "/mtcojo_summary_data.list")
	write.table(gwas, file=gwas_path, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	
	# ref-ld-chr; w-ld-chr
	
	# gwas-thresh
	gsmr_pthres = 5e-8
	
	# mtcojo cmd
	## increase gsmr_pthres until there are more than 10 ind.sig.snps for GSMR
	cond = "1"
	while (cond=="1")
	{
		cmd = paste(gcta_bin, "--threads 10 --diff-freq 0.3 --bfile", ref_geno, "--mtcojo-file", gwas_path, "--ref-ld-chr", ld_chr, "--w-ld-chr", ld_chr, "--out", save_path, "--gwas-thresh", gsmr_pthres)
		system(cmd)
		gsmr_pthres = gsmr_pthres * 100
		cond = system(paste0("grep -v 'nan' ", save_path, ".mtcojo.cma | wc -l"), intern = TRUE)
	}
	
	# sumstats
	res = fread(paste0(save_path, ".mtcojo.cma"), header=FALSE)[, c(1:4,9:11,8)]
	colnames(res) = c("SNP", "EFFECT_ALLELE", "NON_EFFECT_ALLELE", "FRQ", "BETA", "SE", "P", "N")
	res$Z = res$BETA / res$SE
	path_formatted_gwas = Standardise_GWAS(gwas_sumstats_path = res, save_path = paste0(save_path, ".txt.gz"), ref_genome = "GRCh37", convert_ref_genome = NULL)
	
	# sig.mtcojo.snps
	res = fread(path_formatted_gwas)	
	sig.res = res[res$P<pthres, ]
	sig.res = merge(sig.res, merged_path)
	sig.res = sig.res[order(sig.res$P, decreasing=FALSE),]
	fwrite(sig.res, file=sig_save_path, sep="\t", quote=FALSE)

	return(paste0(save_path, ".txt.gz"))
}


