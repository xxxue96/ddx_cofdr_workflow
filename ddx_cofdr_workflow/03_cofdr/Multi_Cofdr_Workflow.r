# run cofdr among traits sharing at least one locus
# run hyprcoloc among all traits (with or without shared locus), and hyprcoloc will cluster colocalized traits through Bayesian divisive clustering algorithm.
# compare cofdr and hyprcoloc results to find out validated colocalized loci
 
source("01_data/00_Process_Sumstats.r")
source("03_cofdr/01_02_Run_Cofdr.r")
source("03_cofdr/03_01_SNP_Split.r")
source("03_cofdr/03_02_Run_Hyprcoloc.r")
source("03_cofdr/04_SNP_Compare.r")
source("03_cofdr/05_01_Run_Gene_Based_Cofdr_MAGMA.r")
source("03_cofdr/05_02_Run_Gene_Based_Cofdr_PrediXcan.r")
source("03_cofdr/05_03_Run_Gene_Based_Cofdr_PWAS.r")
source("03_cofdr/05_04_Gene_Compare.r")
source("03_cofdr/Multi_Cofdr_Workflow_params.r")

#@t2_cofdr: traits to be included in multi-trait cofdr analysis
#@t2_hyprcoloc: traits to be included in hyprcoloc analysis
#@cofdr_name: set a pseudo-name for the multi-cofdr/hyprcoloc derived trait sumstats. Default is NULL to directly concatenate t1 and t2_cofdr as a pseudo-name. 
Run_Multi_Cofdr <- function(t1, t2_cofdr, t2_hyprcoloc, cofdr_name=NULL)
{	
	# basic info for cofdr
	traits = c(t1, t2_cofdr)
	if(is.null(cofdr_name)){ pheno = paste0(traits, collapse="_") } else{ pheno = cofdr_name }
	gwasFilePaths = paste0("01_data/00_Standardise_GWAS/", traits, ".txt.gz")
	merged_path = merge_sumstats(gwasFilePaths, traits)
	# cor.mat = as.matrix(ldsc_intercept_src[traits, traits])
	cor.mat = NULL
	
	# basic info for hyprcoloc
	if(all(t2_hyprcoloc %in% t2_cofdr)){
		traits_hyprcoloc = traits
		gwasFilePaths_hyprcoloc = gwasFilePaths
		merged_path_hyprcoloc = merged_path
	} else{
		traits_hyprcoloc = c(t1, t2_hyprcoloc)
		gwasFilePaths_hyprcoloc = paste0("01_data/00_Standardise_GWAS/", traits_hyprcoloc, ".txt.gz")
		merged_path_hyprcoloc = merge_sumstats(gwasFilePaths_hyprcoloc, traits_hyprcoloc)
	}
		
	# run cofdr
	cofdr_gwas = Run_Cofdr(traits, gwasFilePaths, merged_path, cor.mat=NULL, save_path=paste0("03_cofdr/01_Cofdr_Result/", pheno, ".txt.gz"), sig_save_path=paste0("03_cofdr/02_Cofdr_SigSNP/", pheno, "_", TTthres, "TT.txt.gz"), fuma_save_path=paste0("03_cofdr/05_01_Cofdr_FUMA_Input/", pheno, ".txt.gz"))
	cofdr_snp = cofdr_gwas$sig
			
	# gwas-pw not applicable for multi-trait comparisons
	tmp <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("SNP", "chunkPPA_3", "snpPPA_3"))
	fwrite(tmp, file=paste0(tempdir(), "/gwaspw.tmp.txt"), sep="\t")
	gwas_pw_snp = list(sig=paste0(tempdir(), "/gwaspw.tmp.txt"))
		
	# run hyprcoloc
	tryCatch({
		gwas_split_paths = SNP_Split_preMerge(merged_path_hyprcoloc, traits_hyprcoloc, region_path=ldetect_data, return_data=TRUE, save_path = tempdir())
		hyprcoloc_res = Run_HyPrColoc(traits_hyprcoloc, gwas_list=gwas_split_paths, binary.outcomes, save_path=paste0("03_cofdr/08_HyPrColoc_Result/", pheno, ".txt.gz"))
		# select iterations with analyzed phenotypes containing t2_cofdr	
		tmp = as.data.frame(fread(hyprcoloc_res))
		res=c()
		iter_trait = lapply(tmp$traits, strsplit, ", ")
		for(i in 1:length(iter_trait)){ res =c(res, sum(traits %in% unlist(iter_trait[[i]])) == length(traits))}
		tmp = tmp[res, ] 
		fwrite(tmp, file=paste0(tempdir(), "/hyprcoloc.tmp.txt"), sep="\t")
		hyprcoloc_res = paste0(tempdir(), "/hyprcoloc.tmp.txt")
	}, error=function(e){cat("Run_HyPrColoc Warning :",conditionMessage(e), "\n")})
	
	# compare cofdr and hyprcoloc
	tryCatch({
		sigSNP_comparison = SNP_Compare(snp_path=list(cofdr_snp, gwas_pw_snp$sig, hyprcoloc_res), annot_sumstats=merged_path, save_path=paste0("03_cofdr/09_sigSNP_Comparison/", pheno, ".txt.gz"))
	}, error=function(e){cat("SNP_Compare Warning :",conditionMessage(e), "\n")})

	# basic for gene-based cofdr
	twas_dir = paste0("03_cofdr/04_Cofdr_Gene_Result/", pheno); dir.create(twas_dir, showWarnings = FALSE, recursive = TRUE)
	magma_dir=paste0(twas_dir, "/01_MAGMA"); dir.create(magma_dir, showWarnings = FALSE, recursive = TRUE)
	pwas_dir = paste0(twas_dir, "/04_pwas"); dir.create(pwas_dir, showWarnings = FALSE, recursive = TRUE)

	cofdr_gene_magma = Run_Gene_Based_Cofdr_MAGMA(traits, gwasFilePaths, merged_path, cor.mat, magma_dir, save_dir=twas_dir)
	cofdr_gene_metaxcan = Run_Gene_Based_Cofdr_MetaXcan(traits, gwasFilePaths, merged_path, cor.mat, save_dir=twas_dir)	
	cofdr_gene_pwas = Run_Gene_Based_Cofdr_PWAS(traits, gwasFilePaths, merged_path, cor.mat, pwas_dir, save_dir=twas_dir)
	tryCatch({
		sigGene_comparison = Gene_Compare(traits, magma_path=cofdr_gene_magma, twas_path=cofdr_gene_metaxcan$eqtl, sptwas_path=cofdr_gene_metaxcan$sqtl, pwas_path=cofdr_gene_pwas, save_path=paste0(twas_dir,"/",pheno,"_geneCompare.txt.gz"))
	}, error=function(e){cat("Gene_Compare Warning :",conditionMessage(e), "\n")})
	
	file.remove(paste0(twas_dir,"/01_harmonized_gwas/",traits, ".txt.gz"))
}

# # colocalization analyses among multiple traits with at least one shared locus with COVID-19
# covidhgi = c("A2", "B2", "C2")
# respiratory = c("Asthma", "ILD", "IPF", "COPD", "Pneumonia.meta", "Pneumonia.b.meta", "Pneumonia.v.meta", "Flu.meta")


# TTthres = 0.9
# Run_Multi_Cofdr(t1="A2", t2_cofdr=c("Asthma", "ILD", "IPF", "COPD"), t2_hyprcoloc=respiratory)

# TTthres = 0.8
# Run_Multi_Cofdr(t1="B2", t2_cofdr=c("Asthma", "ILD", "IPF"), t2_hyprcoloc=respiratory)
# Run_Multi_Cofdr(t1="C2", t2_cofdr=c("Asthma", "ILD", "IPF"), t2_hyprcoloc=respiratory)
