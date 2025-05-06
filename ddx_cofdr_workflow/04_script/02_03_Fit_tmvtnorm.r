library(tmvtnorm)
library(data.table)
setDTthreads(40)
source("01_data/00_Process_Sumstats.r")

# MLE fit Truncated Multivariate Normal Distribution(tmvtnorm) to get correlation under the null (i.e. correlation due to sample overlap).
fit.tmvtnorm <- function(gwasFilePaths, traits, merged_path=NULL)
{ 
	common_snp = merge_sumstats(gwasFilePaths, traits, merged_path)
	zmat = as.matrix(common_snp[, c("SNP", grep("^Z_", colnames(common_snp), value=TRUE)), with=FALSE], rownames="SNP")
	zmat = zmat[rowSums(is.na(zmat))==0,]
	
 	# subset zmat with abs(Z)<2 to remove strongly associated snps, and remained snps are very weakly associated so that are treated as snps under the null (i.e. not associated with both phenotypes)
	indices = apply(zmat, 1, function(snp) all(abs(snp) < 2))
	zmat.fit = zmat[indices, ]
	
	# apply mle to fit truncated multivariate normal distribution
	lower = c(-2, -2)
	upper = c(2, 2)
	mle.fit = mle.tmvnorm(zmat.fit, lower=lower, upper=upper)
	
	# extract covariance from fitted model and cal. corrleation 
	# cor.xy = cov.xy / sqrt(cov.xx * cov.yy)
	cor_12 = as.numeric(coef(mle.fit)["sigma_1.2"] / sqrt(coef(mle.fit)["sigma_1.1"] * coef(mle.fit)["sigma_2.2"]))
	logL = as.numeric(logLik(mle.fit))
	
	return(c(cor_12, logL))
}

mynames = commandArgs(TRUE)
mypaths = paste0("01_data/00_Standardise_GWAS/", mynames, ".txt.gz")
overlap_cor = log_Lik = data.frame(matrix(nrow = length(mynames), ncol = length(mynames), dimnames = list(mynames, mynames)))

for(i in mynames){
	for(j in mynames){
		tryCatch(
			if(match(i, mynames) < match(j, mynames)){
				print(paste(i,j, "start"))
				if(is.na(overlap_cor[i,j])){
					traits=c(i,j)
					res = fit.tmvtnorm(gwasFilePaths=as.list(paste0("01_data/00_Standardise_GWAS/", traits, ".txt.gz")), traits)
					overlap_cor[i,j] = res[1]
					log_Lik[i,j] = res[2]
				}
				print(paste(i,j, "done"))
			} else if(match(i, mynames) == match(j, mynames)){
				print(paste(i,j, "start"))
				overlap_cor[i,j] = 1
				log_Lik[i,j] = NA
				print(paste(i,j, "done"))
			} else{
				print(paste(i,j, "start"))
				overlap_cor[i,j] = overlap_cor[j,i]
				log_Lik[i,j] = log_Lik[j,i]
				print(paste(i,j, "done"))
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
}

write.table(t(overlap_cor), file = "01_data/02_tmvtnorm_Results/tmvtnorm_overlapCor.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE, na = "NA")
write.table(t(log_Lik), file = "01_data/02_tmvtnorm_Results/tmvtnorm_overlapCor_logLik.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE, na = "NA")
