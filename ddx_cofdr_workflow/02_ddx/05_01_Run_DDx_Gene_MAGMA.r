library(data.table)
source("04_script/05_MAGMA_Functions.r")
source("02_ddx/DDx_Workflow_params.r")

# step1. run magma 
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

# step2. extract sig.genes
extract_gene_magma <- function(genes_out_path, save_dir)
{
	dat = fread(genes_out_path)[, c("GENE", "P")]
	dat$FDR = p.adjust(dat$P, method="fdr")	
	set.seed(12345)
	dat$Z = abs(qnorm(dat$P/2)) * sample(x=c(-1,1), size=length(dat$P), replace=TRUE)
	
	## filter and annotate sig.genes
	gene_annot_loc = fread(gene_annot)
	dat = merge(dat, gene_annot_loc, by.x="GENE", by.y="entrezID")
	dat = unique(dat[, c("hgnc_symbol", "P", "FDR", "Z")])
	dat = dat[dat$FDR<FDRthres,]
	dat = dat[order(dat$FDR, decreasing=FALSE),]	

	## save
	save_path = paste0(save_dir, "/", gsub(".genes.out", "", basename(genes_out_path)), "_sig.magma.txt.gz")
	fwrite(dat, file=save_path, sep="\t")
	
	return(save_path)
}


Run_DDx_Gene_MAGMA <- function(decor_gwas_paths, magma_dir, save_dir)
{
	magmaPaths = create_magma_gene_results(decor_gwas_paths, magma_dir)
	ddx_gene = extract_gene_magma(magmaPaths[[1]], save_dir)
	
	return(ddx_gene)
}