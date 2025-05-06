# extract h2, rg, ldsc intercept and its pval
library(data.table)

pheno_name = commandArgs(TRUE)
rg = overlap_cor = overlap_cor_pval = data.frame(matrix(nrow = length(pheno_name), ncol = length(pheno_name), dimnames = list(pheno_name, pheno_name)))

for(i in pheno_name){
	for(j in pheno_name){
		if(match(i, pheno_name) < match(j, pheno_name)){
			dat <- tryCatch({as.data.frame(fread(paste("01_data/01_LDSC_Results/rg/",i,"_",j,".log", sep=""), skip="Summary of Genetic Correlation Results"))}, error = function(cond) {message(conditionMessage(cond)); NULL})
			if(!is.null(dat)){
				rg[i,j] = dat$rg
				overlap_cor[i,j] = dat$gcov_int
				overlap_cor_pval[i,j] = 2*pnorm(abs(dat$gcov_int/dat$gcov_int_se), lower.tail=FALSE)
			}			
		} else if(match(i, pheno_name) == match(j, pheno_name)){
			rg[i,j] = NA
			overlap_cor[i,j] = 1
			overlap_cor_pval[i,j] = 0
		} else{
			rg[i,j] = rg[j,i]
			overlap_cor[i,j] = overlap_cor[j,i]
			overlap_cor_pval[i,j] = overlap_cor_pval[j,i]
		}
	}
}

write.table(t(rg), file = "01_data/01_LDSC_Results/ldsc_rg.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE, na = "NA")
write.table(t(overlap_cor), file = "01_data/01_LDSC_Results/ldsc_intercept.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE, na = "NA")
write.table(t(overlap_cor_pval), file = "01_data/01_LDSC_Results/ldsc_intercept_pval.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE, na = "NA")

# negative h2 may cause NA rg https://github.com/bulik/ldsc/issues/88
h2 = data.frame(matrix(nrow = length(pheno_name), ncol = 1, dimnames = list(pheno_name, "h2")))
for(i in pheno_name){
	dat <- as.data.frame(fread(paste("01_data/01_LDSC_Results/h2/",i,"_h2.log", sep=""), skip="Using two-step estimator with cutoff at 30.", header=F, sep=":"))
	h2[i,] = dat[1,2]
}
write.table(h2, file = "01_data/01_LDSC_Results/ldsc_h2.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE, na = "NA")
