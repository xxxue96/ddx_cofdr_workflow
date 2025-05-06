# CCGWAS performs a case-case association testing of two different disorders based on their respective case-control GWAS results https://github.com/wouterpeyrot/CCGWAS
library(CCGWAS)
library(data.table)
setDTthreads(40)
source("02_ddx/DDx_Workflow_params.r")

# @ Neff_param: median(beta_viaNeff/beta_viaOR) reported by CCGWAS error https://github.com/wouterpeyrot/CCGWAS/issues/2
CCGWAS_Input_Reformat <- function(trait, path, Neff, Neff_param=NA)
{
	gwas = fread(path)
	
	if("N_EFF" %in% colnames(gwas)){ colnames(gwas)[which(colnames(gwas)=="N_EFF")] = "Neff" }
	if(!("Neff" %in% colnames(gwas))){ gwas$Neff = rep(as.numeric(Neff), nrow(gwas)) }
	
	if(!is.na(Neff_param)){
		print(paste(basename(path), "adjust Neff with", Neff_param))
		gwas$Neff = Neff_param^2 * gwas$Neff
	} else{
		print(paste(basename(path), "Neff with no adjustment"))
	}
	
	if(! "EA" %in% colnames(gwas)) {gwas$EA = gwas$A2}
	if(! "NEA" %in% colnames(gwas)) {gwas$NEA = gwas$A1}
	if(! "OR" %in% colnames(gwas)) {gwas$OR = exp(gwas$BETA)}
	gwas= gwas[, c("SNP", "CHR", "BP", "EA", "NEA", "FRQ", "OR", "SE", "P", "Neff")]
	
	formatted_path = paste0(tempdir(), "/", basename(path))
	fwrite(gwas, file=formatted_path, sep="\t")
	return(formatted_path)
}

# @save_path: output file path of ccgwas results
# @prevalence_path: file path of phenotype prevenlence.
# @size_path: file path of sample size N, number of cases N_CAS and number of controls N_CON
# @h2_path: path of heritability on liability scale. 
# @rg_path: path of genetic correlation
# @intercept_path: path of ldsc intercept
# @m: approximation of number of independent effective loci
# @subtype_data: FALSE when comparing two different disorders. TRUE when compating subtypes of a disorder
CCGWAS_Function <- function(t1, t2, path1, path2, save_path, 
							prevalence_path=NULL, 
							size_path=sample_size_path,
							h2_path=ldsc_obs_h2_path,
							rg_path=ldsc_rg_path,
							intercept_path=ldsc_intercept_path,
							m=5000, subtype_data=FALSE)
{
	## sample size info
	size = read.table(size_path, row.names=1, header=TRUE, sep="\t")
	N_A1 = size[t1, "N_CAS"]; N_A0 = size[t1, "N_CON"]; N_A =  size[t1, "N"]
	N_B1 = size[t2, "N_CAS"]; N_B0 = size[t2, "N_CON"]; N_B =  size[t2, "N"]
	N_overlap_A0B0 = 0
	
	## prevenlence info
	## prevenlence=N_CAS/N, confidence interval for the propertion was then estimated using Wald interval https://stats.stackexchange.com/questions/568892/calculate-a-95-confidence-interval-of-a-population-proportion-in-r.
	if(is.null(prevalence_path)){
		K_A1A0 = N_A1/N_A; K_B1B0 = N_B1/N_B
		CI.A = K_A1A0 + qnorm(c(.025,.975))*sqrt(K_A1A0*(1-K_A1A0)/N_A); K_A1A0_low = CI.A[1]; K_A1A0_high = CI.A[2]
		CI.B = K_B1B0 + qnorm(c(.025,.975))*sqrt(K_B1B0*(1-K_B1B0)/N_B); K_B1B0_low = CI.B[1]; K_B1B0_high = CI.B[2]
	} else{
		prevenlence = read.table(prevalence_path, row.names=1, header=TRUE, sep="\t")
		K_A1A0 = prevalence[t1, "K"]; K_A1A0_low = prevalence[t1, "K_low"]; K_A1A0_high = prevalence[t1, "K_high"]
		K_B1B0 = prevalence[t2, "K"]; K_B1B0_low = prevalence[t2, "K_low"]; K_B1B0_high = prevalence[t2, "K_high"]
	}
	
	## ldsc info
	## assume sample prevenlence = population prevenlence, and cal. liablity heritability.
	ldsc_h2 = read.table(h2_path, row.names=1, header=TRUE, sep="\t")
	ldsc_rg = read.table(rg_path, row.names=1, header=TRUE, sep="\t")
	ldsc_intercept = read.table(intercept_path, row.names=1, header=TRUE, sep="\t")
	h2l_A1A0 = as.numeric(strsplit(ldsc_h2[t1, "h2"], " ")[[1]][1])
	h2l_B1B0 = as.numeric(strsplit(ldsc_h2[t2, "h2"], " ")[[1]][1])
	rg_A1A0_B1B0 = ldsc_rg[t1,t2]
	intercept_A1A0_B1B0 = ldsc_intercept[t1, t2]
	
	## sumstats reformat
	Neff_A1A0 = size[t1, "N_EFF"]; Neff_A1A0_param = size[t1, "N_EFF_param"]
	Neff_B1B0 = size[t2, "N_EFF"]; Neff_B1B0_param = size[t2, "N_EFF_param"]
	path1 = CCGWAS_Input_Reformat(t1, path1, Neff_A1A0, Neff_A1A0_param)
	path2 = CCGWAS_Input_Reformat(t2, path2, Neff_B1B0, Neff_B1B0_param)
	
	## run ccgwas
	CCGWAS(outcome_file=save_path, A_name=t1, B_name=t2, sumstats_fileA1A0=path1, sumstats_fileB1B0=path2,
		   K_A1A0, K_A1A0_high, K_A1A0_low, K_B1B0, K_B1B0_high, K_B1B0_low,
		   h2l_A1A0, h2l_B1B0, rg_A1A0_B1B0, intercept_A1A0_B1B0, m,
		   N_A1, N_B1, N_A0, N_B0, N_overlap_A0B0, subtype_data)
		
	save_path = paste0(save_path, ".results.gz")
	return(save_path)
}

