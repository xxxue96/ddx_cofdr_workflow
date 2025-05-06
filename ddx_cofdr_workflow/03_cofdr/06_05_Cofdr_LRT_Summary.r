# summarize prop and LRT for each comparison
library(data.table)

Cofdr_LRT_Summary <- function(t1, t2)
{
	pheno = paste0(t1, "_", t2)

	snp = fread(paste0("03_cofdr/01_Cofdr_Result/", pheno, ".log"))[, c("piTT", "LRTp")]	
	magma = fread(paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno, "/01_MAGMA/", pheno, ".0UP.0DOWN.log"))[, c("piTT", "LRTp")]	
	twas_lung = fread(paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno, "/02_spredixcan/eqtl/", pheno, "_Lung.log"))[, c("piTT", "LRTp")]	
	twas_whole_blood = fread(paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno, "/02_spredixcan/eqtl/", pheno, "_Whole_Blood.log"))[, c("piTT", "LRTp")]	
	twas_smultixcan = fread(paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno, "/03_smultixcan/eqtl/", pheno, "_smultixcan.log"))[, c("piTT", "LRTp")]	
	sptwas_lung = fread(paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno, "/02_spredixcan/sqtl/", pheno, "_Lung.log"))[, c("piTT", "LRTp")]	
	sptwas_whole_blood = fread(paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno, "/02_spredixcan/sqtl/", pheno, "_Whole_Blood.log"))[, c("piTT", "LRTp")]	
	sptwas_smultixcan = fread(paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno, "/03_smultixcan/sqtl/", pheno, "_smultixcan.log"))[, c("piTT", "LRTp")]	
	pwas = fread(paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno, "/04_pwas/", pheno, ".log"))[, c("piTT", "LRTp")]
	
	res = cbind(pheno, snp, magma, twas_lung$piTT, twas_whole_blood$piTT, twas_smultixcan$piTT, twas_lung$LRTp, twas_whole_blood$LRTp, twas_smultixcan$LRTp, sptwas_lung$piTT, sptwas_whole_blood$piTT, sptwas_smultixcan$piTT, sptwas_lung$LRTp, sptwas_whole_blood$LRTp, sptwas_smultixcan$LRTp, pwas)
	colnames(res) = c("Comparisons", "prop.SNP", "LRT.SNP", "prop.Gene.MAGMA", "LRT.Gene.MAGMA", paste0("prop.Gene.TWAS.", c("lung", "whole_blood", "smultixcan")), paste0("LRT.Gene.TWAS.", c("lung", "whole_blood", "smultixcan")), paste0("prop.Gene.spTWAS.", c("lung", "whole_blood", "smultixcan")), paste0("LRT.Gene.spTWAS.", c("lung", "whole_blood", "smultixcan")), "prop.Gene.PWAS", "LRT.Gene.PWAS")
	
	return(res)
}

# covidhgi = c("A2", "B2", "C2")
# respiratory = c("Asthma", "ILD", "IPF", "COPD", "Pneumonia.meta", "Pneumonia.b.meta", "Pneumonia.v.meta", "Flu.meta")

# res = data.frame()
# for(t1 in covidhgi){
	# res = rbind(res, do.call(rbind.data.frame, lapply(respiratory, Cofdr_LRT_Summary, t1=t1)))
# }

# fwrite(res, file="03_cofdr/10_Excel_Summary/Tables/COVID_Cofdr_LRT_Summary.txt", sep="\t")