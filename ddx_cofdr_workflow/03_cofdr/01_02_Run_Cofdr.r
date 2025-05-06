# @traits: pheontypes of analyzed gwas
# @gwasFilePaths: gwas sumstats for each phenotype to run colocalization. paths are also accepted.
# @merged_path: merged gwas sumstats. paths are also accepted.
# @cor.mat: sample overlap correlation matrix.
# @save_path: full path to save cofdr results

library(data.table)
setDTthreads(20)
source("01_data/00_Process_Sumstats.r")
source("03_cofdr/01_00_Zmat_Decor.r")
source("03_cofdr/01_01_Cofdr_Functions.r")
source("03_cofdr/Cofdr_Workflow_params.r")

Run_Cofdr <- function(traits, gwasFilePaths, merged_path=NULL, cor.mat=NULL, save_path, sig_save_path=NULL, fuma_save_path=NULL, keep_Zstat=FALSE)
{	
	##  decorrelate zmat
	common_snp = merge_sumstats(gwasFilePaths, traits, merged_path)
	zmat = as.matrix(common_snp[, c("SNP", grep("^Z_", colnames(common_snp), value=TRUE)), with=FALSE], rownames="SNP")
	if(!is.null(cor.mat)) {zmat = t(decor(t(zmat), cor.mat))}
	
	## run cofdr	
	result = cofdr(N=ncol(zmat), zmat, iter=500, tol=1e-8, empNULL=FALSE, nulltype=1, bre=120, df=7)
	post.prob.mat = cbind(result$post.prob[, c(1:3,ncol(result$post.prob))], common_snp) # snp significant in at most 1 trait, or in all traits	
	colnames(post.prob.mat) = c("FF","TF","FT","TT",colnames(common_snp))
	post.prob.mat = post.prob.mat[order(post.prob.mat$TT, decreasing=TRUE), ]
	if(keep_Zstat==TRUE){
		fwrite(post.prob.mat, file=save_path, sep="\t", quote=FALSE)
	} else{
		fwrite(post.prob.mat[, c("FF","TF","FT","TT","SNP")], file=save_path, sep="\t", quote=FALSE)
	}
	
	## cofdr stats: pi.vector, diff, D, LRT.p
	stat = as.data.table(t(c(result$pi.vector[c(1:3,ncol(result$post.prob))], result$diff, result$D, result$LRT.p)))
	stat = cbind(paste(traits, collapse = "_"), stat)
	colnames(stat) = c("trait", "piFF", "piTF", "piFT", "piTT", "diff", "D", "LRTp")	
	fwrite(stat, file=gsub(".txt.gz", ".log", save_path), sep="\t", quote=FALSE)
	
	## extract sig.snp
	post.prob.mat$P = 1 - post.prob.mat$TT
	sig.snp = post.prob.mat[post.prob.mat$TT >= TTthres, ]
	if(!is.null(sig_save_path)){fwrite(sig.snp, file=sig_save_path, sep="\t", quote=FALSE)}
	
	## prepare for FUMA
	if(!is.null(fuma_save_path) & nrow(sig.snp)>0){
		fuma_scale = 1e-5/(1-TTthres)
		post.prob.mat$P_FUMA = (1 - post.prob.mat$TT) * fuma_scale
		if(!("N" %in% colnames(post.prob.mat))) {post.prob.mat$N = rowSums(post.prob.mat[, paste0("N_", traits), with=FALSE])}
		fwrite(post.prob.mat[, c("SNP","A2","TT","P_FUMA","N")], file=fuma_save_path, sep="\t", quote=FALSE, scipen=999)
	}
	
	return(list(allsnp=save_path, sig=sig_save_path, fuma=fuma_save_path))
}