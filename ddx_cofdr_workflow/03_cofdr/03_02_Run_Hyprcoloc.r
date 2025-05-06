# https://rdrr.io/github/jrs95/hyprcoloc/f/vignettes/hyprcoloc.Rmd
# Assumption/default parameters: 1) all studies have the same participants, i.e. sample.overlap = matrix(rep(1, dim(effect.est)[2]^2)
#			 					 2) independence between the studies and the traits, i.e. trait.cor = diag(1, dim(effect.est)[2]
#								 3) indenpendence between snps, i.e. ld.matrix = diag(1, dim(effect.est)[1]
#								 4) single causal variant <- Like coloc, hyprcoloc use fine-mapping to find out a candidate snp to explain shared association
# When analysing correlated traits, assuming independence between the studies not only correctly identifies the three clusters of colocalized traits, but does so around 3500 times faster than accounting for non-indepdence.
# i.e. using default settings of sample.overlap, trait.cor, and ls.matrix will give us a correct estimate when analysing correlated traits.

# posterior_prob: the posterior probability that these traits are colocalized
# Bayesian divisive clustering algorithm: when no shared variants among all traits, use the algorithm automatically subset traits to most confident colocalized clusters  	
	# regional_prob: P(R), PPA for subset traits as colocalized clusters using regional selection criterion
	#			  -> computed from a collection of hypotheses which assume that all traits do not colocalize because one of the traits does not have a causal variant in the region
	# align_prob: P(A),  PPA for subset traits as colocalized clusters using alignment selection criterion
	#			  ->  computed from hypotheses which assume that all traits do not colocalize because one of the traits has a causal variant elsewhere in the region
	# considering computing stress and time consuming, regional_prob is default subsetting threthold	
# candidate_snp: a candidate causal variant explaining the shared association -> single causal variant assumption
# posterior_explained_by_snp: the proportion of the posterior probability explained by this variant (which represents the HyPrColoc multi-trait fine-mapping probability)

# P(RA) vs P(R): posterior_prob v.s. regional_prob
#				 P(RA)=P(R)*P(A)<P(R); 
#				 P(R) is used to determine whether traits will be subsetted into one cluster while P(RA) is used to determine whether clustered traits are indeed colocalized
#				 -> outputted P(R) is the most approriate value to get maximum P(RA) -> sometimes P(R)>0.9 while P(RA)=0 may appear, which means that clustered traits are not colocalized anyway.

library(tidyverse)
library(data.table)
library(hyprcoloc)
setDTthreads(20)	
source("01_data/00_Process_Sumstats.r")

# @traits: trait{1...N}
# @gwas_list: gwas list for trait{1...N}. Each gwas with columns "SNP, BETA_trait{1...N}, SE_trait{1...N}, chunk" correspondingly
Run_HyPrColoc <- function(traits, gwas_list, binary.outcomes, save_path)
{
	if(length(gwas_list)==1 | is.data.frame(gwas_list)){
		gwas_list = read_sumstats(gwas_list)
	} else{
		gwas_list = lapply(gwas_list, read_sumstats)
		gwas_list = gwas_list %>% reduce(inner_join, by = c("SNP", "chunk"))
	}
	
	run_hyprcoloc_per_region <- function(dat)
	{
		tryCatch({
			rsid = dat$SNP
			
			betas = data.matrix(dat[, paste0("BETA_", traits)])
			rownames(betas) = rsid
			colnames(betas) = traits

			ses = data.matrix(dat[, paste0("SE_", traits)])
			rownames(ses) = rsid
			colnames(ses) = traits
			
			res = hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid, binary.outcomes=binary.outcomes, snpscores=FALSE)
			
			# filter colocalized results
			res = res$results[res$results$traits != "None", ]
			return(res)
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}
	gwas_list = split(gwas_list, f = gwas_list$chunk)
	hyprcoloc.res = do.call(rbind.data.frame, lapply(gwas_list, run_hyprcoloc_per_region))
	
	fwrite(hyprcoloc.res, file = save_path, sep = "\t", quote=FALSE)
	return(save_path)
}


























