# Split SNPs into LD-independent genomic regions
library(data.table)
setDTthreads(20)

# @path: gwas sumstats path
# @region_path: ld-independent path
# @save_path: path to save split gwas results
SNP_Split <- function(path, region_path, return_data = TRUE, save_path = tempdir())
{
	gwas = fread(path)[, c("CHR", "BP", "BP", "SNP", "BETA", "SE", "N")]
	gwas_header = colnames(gwas)
		
	region = fread(region_path)
	region$chunk = as.numeric(rownames(region))
	region$chr = gsub(pattern="chr", "", region$chr)
	region_header = colnames(region)

	library(bedtoolsr)
	gwas_sort = bt.sort(i = gwas)
	region_sort = bt.sort(i = region)
	gwas_split = bt.intersect(a = gwas_sort, b = region_sort, wo = "")	# SNP locating exactly at the chunk_start/stop will be counted twice e.g. rs7546735
	
	gwas_split = gwas_split[, c(4:6,11)]
	pheno = gsub(".txt.gz", "", basename(path))
	colnames(gwas_split) = c("SNP", paste0("BETA_",pheno), paste0("SE_",pheno), "chunk")
	
	if(return_data == TRUE){
		return(gwas_split)
	} else{
		save_path = paste0(save_path, "/", basename(path))
		fwrite(gwas_split, file = save_path, sep = " ")
		return(save_path)
	}
}

# directly split paths of pre-merged SNPs
# @merged_path: merged sumstats of multiple traits
SNP_Split_preMerge <- function(merged_path, traits, region_path, return_data = TRUE, save_path = tempdir())
{
	gwas = read_sumstats(merged_path)[, c("CHR", "BP", "BP", "SNP", paste0("BETA_",traits), paste0("SE_",traits)), with=FALSE]
	gwas_header = colnames(gwas)
		
	region = fread(region_path)
	region$chunk = as.numeric(rownames(region))
	region$chr = gsub(pattern="chr", "", region$chr)
	region_header = colnames(region)
	
	# old_path = Sys.getenv("PATH")
	# Sys.setenv(PATH = paste(old_path, "/home/xue/Python/anaconda3/pkgs/bedtools-2.30.0-hc088bd4_0/bin", sep=":"))
	library(bedtoolsr)
	gwas_sort = bt.sort(i = gwas)
	region_sort = bt.sort(i = region)
	gwas_split = bt.intersect(a = gwas_sort, b = region_sort, wo = "")	# SNP locating exactly at the chunk_start/stop will be counted twice e.g. rs7546735
	colnames(gwas_split) = c(colnames(gwas), colnames(region))
	
	gwas_split = gwas_split[,  c("SNP", paste0("BETA_",traits), paste0("SE_",traits), "chunk")]
	
	if(return_data == TRUE){
		return(gwas_split)
	} else{
		save_path = paste0(save_path, "/", traits)
		fwrite(gwas_split, file = save_path, sep = " ")
		return(save_path)
	}
}