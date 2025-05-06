library(data.table)
setDTthreads(20)
source("01_data/00_Process_Sumstats.r")
source("02_ddx/01_sddx_Function.r")
source("02_ddx/DDx_Workflow_params.r")

# @ldsc_intercept: LD score regression intercept, numeric.
Run_DDx <- function(t1, t2, path1, path2, merged_path=NULL, ldsc_intercept, save_path, sig_save_path)
{	
	# merge two gwas and extract common snps with the same order
	# merged_path: whether merged results are avaliable. Default is NA to merge gwas within the function
	common_snp = merge_sumstats(gwasFilePaths=list(path1,path2), traits=c(t1,t2), merged_path)
	
	gwas1 = common_snp[, c("SNP","CHR","BP","A1","A2", paste0(c("FRQ","BETA","SE","P","Z","N"), "_", t1)), with=FALSE]
	gwas2 = common_snp[, c("SNP","CHR","BP","A1","A2", paste0(c("FRQ","BETA","SE","P","Z","N"), "_", t2)), with=FALSE]
	colnames(gwas1) = colnames(gwas2) = c("SNP","CHR","BP","A1","A2","FRQ","BETA","SE","P","Z","N")
	
	# run ddx
	print(paste(t1, t2, "start with overlap correlation", ldsc_intercept))
	ddx.gwas = sddx(ldsc_intercept, gwas1, t1.case=gwas1$N, gwas2, t2.case=gwas2$N)	
	ddx.gwas = ddx.gwas[order(ddx.gwas$P, decreasing=FALSE),]
	fwrite(ddx.gwas, file=save_path, sep="\t", na="NA", quote=FALSE, scipen=999)
	
	# extract sig.ddx.snp
	sig.snp = ddx.gwas[ddx.gwas$P<=pthres, ]
	sig.snp = merge(sig.snp, common_snp, by=c("SNP", "CHR", "BP", "A1", "A2"))
	fwrite(sig.snp, file=sig_save_path, sep="\t", na="NA", quote=FALSE)
	
	return(save_path)
}


