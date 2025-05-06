# clump SNPs based on ld matrix generated from 1000G reference genome using ieugwas api: https://github.com/MRCIEU/TwoSampleMR/blob/d498fb173797aa5fa4b813cf4fde8418e35c601c/R/ld.R
# only bi-allelic SNPs with MAF > 0.01 are included "https://gwaslab.com/2021/04/23/%e8%a7%a3%e9%87%8a%e5%a4%8d%e6%9d%82%e7%96%be%e7%97%85%e7%9a%84%e7%a7%8d%e4%b8%bb%e6%b5%81%e6%a8%a1%e5%9e%8b-cdcv-rame-infinitesimal-broad-sense-heritability/"
library(ieugwasr)
library(data.table)
setDTthreads(40)
source("01_data/00_Process_Sumstats.r")

SNP_Clump <- function(path, save_path, clump_kb=1000, clump_r2=0.01, clump_p=5e-8, pop="EUR", bfile=NULL, plink_bin=NULL, return_data=FALSE)
{	
	dat = read_sumstats(path)
	colnames(dat)[which(colnames(dat)=="SNP")] = "rsid"
	colnames(dat)[which(colnames(dat)=="P")] = "pval"
	
	dat = ieugwasr::ld_clump(dat, clump_kb, clump_r2, clump_p, pop, opengwas_jwt = get_opengwas_jwt(), bfile, plink_bin)	
	dat = dat[, -ncol(dat), with=FALSE]
	colnames(dat)[which(colnames(dat)=="rsid")] = "SNP"
	colnames(dat)[which(colnames(dat)=="pval")] = "P"
	
	if(return_data == TRUE){
		return(dat)
	} else{
		dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
		fwrite(dat, file=save_path, sep="\t")
		return(save_path)	
	}
}
