# Colocalized genes based on magma gene pvalues (map snps to genes based on physical distances)
library(dplyr)
library(data.table)
setDTthreads(20)
source("01_data/00_Process_Sumstats.r")
source("03_cofdr/01_00_Zmat_Decor.r")
source("03_cofdr/01_02_Run_Cofdr.r")
source("04_script/05_MAGMA_Functions.r")
source("03_cofdr/Cofdr_Workflow_params.r")

# step1: decorrelate gwas summary statistics
# @decor_path: save path of decorrelated gwas 
decor_snp4magma <- function(traits, gwasFilePaths, merged_path=NULL, cor.mat=NULL, decor_path=tempdir())
{	
	##  decorrelate zmat and connvert to 2-sided pmat
	common_snp = merge_sumstats(gwasFilePaths, traits, merged_path)
	zmat = as.matrix(common_snp[, c("SNP", grep("^Z_", colnames(common_snp), value=TRUE)), with=FALSE], rownames="SNP")
	if(!is.null(cor.mat)) {zmat = t(decor(t(zmat), cor.mat))}
	pmat = 2*pnorm(-abs(zmat)); colnames(pmat) = paste0("P_", traits)
	
	## save decorrelated gwas with SNP, P, N	
	pmat = cbind(common_snp[, c("SNP", paste0("N_", traits)), with=FALSE], pmat)
	create_decor_gwas_paths <- function(pheno)
	{
		decor_gwas = pmat[, c("SNP", paste0(c("P_", "N_"), pheno)), with=FALSE]
		colnames(decor_gwas) <- c("SNP", "P", "N")
		
		outPath = paste0(decor_path, "/", pheno, ".txt.gz")
		fwrite(decor_gwas, file=outPath, sep="\t")
		return(outPath)
	}
	decor_gwas_paths = lapply(traits, create_decor_gwas_paths)
	
	return(decor_gwas_paths)
}

# step2: run magma
create_magma_gene_results <- function(decor_gwas_paths, magma_dir)
{
	run_magma <- function(gwas_path)
	{
		magmaPaths = get_magma_paths(path=gwas_path, output_path=magma_dir, upstream_kb, downstream_kb, id_type)
		genes_out_path = map_snps_to_genes(path=gwas_path, genome_build="GRCh37", upstream_kb, downstream_kb, N=NULL, genome_ref_path, population=tolower(pop), genes_only=FALSE, storage_dir=tools::R_user_dir("MAGMA.Celltyping", which="cache"), force_new=FALSE, verbose=TRUE,
										   magmaPaths, genomeLocFile, genesAnnotPrefix)
		return(genes_out_path)
	}
	magmaPaths = lapply(decor_gwas_paths, run_magma)
		
	return(magmaPaths)
}

# step3: run gene-based cofdr with zmat obtained from step2
# @magma_dir: directory to save magma results
# @save_dir: directory to save significant cofdr genes
run_magma_cofdr <- function(traits, magmaPaths, magma_dir, save_dir)
{
	## Different from ZSTAT from MAGMA is calculated by qnorm(p, lower.tail=FALSE), we cal. absolute z-score from magma gene-based p-value (converted from a F-test) and assign random +/- (extreme +/- zscore will give the same p-value=pnorm(-abs(z-score),lower.tail=TRUE))
	convert_magma_p_to_z <- function(genes_out_path)
	{
		dat = fread(genes_out_path)
		pheno = gsub(paste0(".", upstream_kb, "UP.", downstream_kb, "DOWN.genes.out"), "", basename(genes_out_path))		
				
		set.seed(nchar(pheno)*nchar(genes_out_path))
		dat$Z = abs(qnorm(dat$P/2)) * sample(x=c(-1,1), size=length(dat$P), replace=TRUE)		
		
		dat = dat[, c("GENE", "Z")]
		colnames(dat) = c("SNP", paste0("Z_",pheno))
		
		return(dat)
	}
	magmaGWAS = lapply(magmaPaths, convert_magma_p_to_z)
	
	## run cofdr for magam gene zscores, remove NA/Inf for zscore cols
	merged_path = magmaGWAS %>% reduce(inner_join, by = "SNP")
	merged_path = merged_path[is.finite(rowSums(merged_path[, -"SNP"])), ]
	pheno = paste(traits, collapse="_")
	save_path = paste0(magma_dir, "/", pheno, ".", upstream_kb, "UP.", downstream_kb, "DOWN.txt.gz")
	gene_based_cofdr = Run_Cofdr(traits, gwasFilePaths=NULL, merged_path, cor.mat=NULL, save_path, keep_Zstat=TRUE)
	gene_based_cofdr = gene_based_cofdr$allsnp
	
	## filter and annotate cofdr results
	gene_based_cofdr = fread(gene_based_cofdr)[, c("SNP", "TT")]
	colnames(gene_based_cofdr)[1] <- "entrezID"
	gene_annot_loc = fread(gene_annot)
	gene_based_cofdr = merge(gene_based_cofdr, gene_annot_loc, by="entrezID")
	gene_based_cofdr = gene_based_cofdr[order(gene_based_cofdr$TT, decreasing=TRUE),]		
	fwrite(gene_based_cofdr, file=save_path, sep="\t")
	
	sig_gene_based_cofdr = gene_based_cofdr[gene_based_cofdr$TT>=TTthres_magma,]
	sig_save_path = paste0(save_dir, "/", pheno, ".", upstream_kb, "UP.", downstream_kb, "DOWN_sig.magma.txt.gz")
	fwrite(sig_gene_based_cofdr, file=sig_save_path, sep="\t")
	
	return(sig_save_path)
}

# wrap-up function
Run_Gene_Based_Cofdr_MAGMA <- function(traits, gwasFilePaths, merged_path, cor.mat=NULL, magma_dir, save_dir)
{
	decor_gwas_paths = decor_snp4magma(traits, gwasFilePaths, merged_path, cor.mat, decor_path=tempdir())
	magmaPaths = create_magma_gene_results(decor_gwas_paths, magma_dir)
	gene_based_cofdr = run_magma_cofdr(traits, magmaPaths, magma_dir, save_dir)
	
	return(gene_based_cofdr)
}
