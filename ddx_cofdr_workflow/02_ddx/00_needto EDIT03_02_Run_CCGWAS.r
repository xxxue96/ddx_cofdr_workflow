library(dplyr)
source("02_ddx/03_01_CCGWAS_Function.r")

# subtype_data=TRUE when comparing subtypes of a disorder (with same definition of controls)
Run_CCGWAS <- function(t1, t2, path1, path2, save_path, subtype_data=FALSE, sig_save_path)
{
	if(is.na(ldsc_rg_src[t1, t2]) | ldsc_rg_src[t1, t2]>0.8){
		sink(paste0(save_path, ".log"))
		cat("CC-GWAS is intended for comparing two different disorders with genetic correlation <=0.8")
		cat("\n")
		cat(paste("rg between", t1, "and", t2, "is", ldsc_rg_src[t1, t2]))
		sink()
	} else{
		error = tryCatch({CCGWAS_Function(t1, t2, path1, path2, save_path, subtype_data=subtype_data)}, error=function(e){return(conditionMessage(e))})
		
		## if there is any error, adjust Neff to bypass double-check and re-run ccgwas https://github.com/wouterpeyrot/CCGWAS/issues/2
		if(error != paste0(save_path, ".results.gz")){
			if(grepl(">", error)){
				sample_size_src[t2, "N_EFF_param"] = as.numeric(qdapRegex::ex_between(error, "=", ">")[[1]])
			} else{
				sample_size_src[t2, "N_EFF_param"] = as.numeric(qdapRegex::ex_between(error, "=", "<")[[1]])
			}
			write.table(sample_size_src, file=sample_size_path, row.names = TRUE, col.names = TRUE, sep = "\t", quote=FALSE)
			CCGWAS_Function(t1, t2, path1, path2, save_path, subtype_data=subtype_data)
		}

		dat = fread(paste0(save_path, ".results.gz"))
		dat = dat %>% rename("A2" = "EA", "A1" = "NEA")
		sig.dat = dat[dat$CCGWAS_signif==1, ]
		sig.dat = sig.dat[order(sig.dat$OLS_pval, decreasing=FALSE), ]
		fwrite(sig.dat, file=sig_save_path, sep="\t", na="NA", quote=FALSE)
	}
	
	return(paste0(save_path, ".results.gz"))
}

