library(data.table)
library(dplyr)
setDTthreads(40)

# run pwas and extract sig.genes
Run_DDx_Gene_PWAS <- function(traits, decor_gwas_paths, pwas_dir, save_dir)
{
	run_pwas <- function(i)
	{
		pheno=traits[i]; INPUT=paste0(tempdir(),"/",pheno,".txt.gz"); OUTPUT=pwas_dir
		# change "A1" as effect allele
		tmp = fread(decor_gwas_paths[[i]])
		colnames(tmp)[which(colnames(tmp) %in% c("A1", "A2"))] = c("A2", "A1")
		fwrite(tmp, file=INPUT, sep="\t")
		cmd = paste("bash 04_script/07_Run_PWAS.bash", pheno, INPUT, OUTPUT)
		system(cmd)
		file.remove(INPUT)
		
		res = paste0(pwas_dir, "/", pheno, "_chr", 1:22, ".out")
		res = do.call(rbind.data.frame, lapply(res, fread))
		res = res[, c("ID", "PWAS.P", "PWAS.Z", "FILE")]
		res$FDR = p.adjust(res$`PWAS.P`, method="fdr")
		res = res[res$FDR<FDRthres,]
		res = res[order(res$FDR, decreasing=FALSE),]
		
		return(res)
	}
	pwasRes = lapply(1:length(traits), run_pwas)
	sig_ddx_pqtl = do.call(rbind.data.frame, pwasRes)
	colnames(sig_ddx_pqtl)[1] = "hgnc_symbol"
	
	path_sig_ddx_pqtl = paste0(save_dir, "/", traits, "_sig.pqtl.txt.gz")
	fwrite(sig_ddx_pqtl[, c("hgnc_symbol", "PWAS.P", "FDR", "PWAS.Z", "FILE")], file=path_sig_ddx_pqtl, sep="\t")
	return(path_sig_ddx_pqtl)
}
