# single causal variant assumption for gwas-pw
library(data.table)
setDTthreads(20)
source("03_cofdr/Cofdr_Workflow_params.r")

# @pheno: phenotype name
# @path: path for corresponding gwas summary statistics 
Read_GWAS <- function(pheno, path)
{
	gwas = fread(path)
	
	# Keep the first duplicated record only
	gwas = gwas[!duplicated(gwas$SNP) & is.finite(gwas$BETA) & is.finite(gwas$SE), ]
	
	# Cal. Zscore and Variance; Change CHR prefix
	CHR_chr <- character(nrow(gwas))
	Z <- V <- numeric(nrow(gwas))
	attach(gwas)
	for (i in 1:nrow(gwas))
	{
		CHR_chr[i] = paste0("chr", CHR[i])
		Z[i] = BETA[i] / SE[i]
		V[i] = SE[i] * SE[i]
	}
	detach(gwas)
	gwas = cbind(gwas, CHR_chr, Z, V)
	gwas = gwas[, c("SNP", "CHR_chr", "BP", "Z", "V")]
		
	# Change column names
	colnames(gwas) = c("SNPID", "CHR", "POS", paste0("Z_",pheno), paste0("V_",pheno))
	
	return(gwas)
}

Run_gwas_pw <- function(t1, t2, path1, path2, save_path, ldetect_data, overlap_cor)
{	
	# check gwas-pw instalation
	# old_path = Sys.getenv("PATH")
	# Sys.setenv(PATH = paste(old_path, "/mnt/data/xue/Tool/gwas-pw-0.21/src", sep=":"))

	# data reformat
	gwas1 = Read_GWAS(t1, path1)
	gwas2 = Read_GWAS(t2, path2)
	
	# merged sumstats should be ordered by chromosome position and no scientific notation is allowed
	common_gwas = merge(gwas1, gwas2, by = c("SNPID", "CHR", "POS"))
	common_gwas$CHR = factor(x = common_gwas$CHR, levels = paste0("chr",1:22), ordered=TRUE)
	common_gwas = common_gwas[order(common_gwas$CHR, common_gwas$POS), ]	
	common_gwas_path = paste0(tempdir(), "/", t1, "_", t2, ".txt.gz")	
	fwrite(common_gwas, file = common_gwas_path, sep = " ", scipen=999)	
		
	# run gwas-pw
	## remove snps outside regions defined in ldetect_data until there is no error	
	gwas_pw_cmd = paste("gwas-pw -i", common_gwas_path, "-bed", ldetect_data, "-phenos", t1, t2, "-o", save_path, "-cor", overlap_cor, "2>&1")
	error = system(gwas_pw_cmd, intern=TRUE)
	print(error)
	
	while (sum(grepl("ERROR", error)) == 1) {
		pos = error[grepl("ERROR", error)]
		pos = as.numeric(sub(".*position is ", "", pos))
		print("remove SNP")
		print(common_gwas[which(common_gwas$POS==pos), ])
		common_gwas = common_gwas[-which(common_gwas$POS==pos), ]
		fwrite(common_gwas, file = common_gwas_path, sep = " ", scipen=999)
		
		gwas_pw_cmd = paste("gwas-pw -i", common_gwas_path, "-bed", ldetect_data, "-phenos", t1, t2, "-o", save_path, "-cor", overlap_cor, "2>&1")
		error = system(gwas_pw_cmd, intern=TRUE)		
	}
	
	return(save_path)
}


