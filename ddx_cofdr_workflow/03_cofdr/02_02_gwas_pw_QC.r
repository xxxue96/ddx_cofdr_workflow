# find out the snp with maximum PPA3 for each chunk
# filter 1447 chunks with chunkPPA3=0.8 to get colocalized regions
library(dplyr)
library(data.table)
setDTthreads(20)

gwas_pw_QC <- function(gwas_pw_res, chunkPPA3, snpPPA3, save_path)
{	
	## find the SNP with max PPA_3 for each chunk
	snpFile = fread(paste0(gwas_pw_res, ".bfs.gz"))
	max_snpPPA = as.data.frame(snpFile %>% group_by(chunk) %>% top_n(1, PPA_3))[, c("id", "PPA_3", "chunk")]
	colnames(max_snpPPA) = c("SNP", "snpPPA_3", "chunk")
	
	## match regional PPA3 for each chunk 
	chunkFile = fread(paste0(gwas_pw_res, ".segbfs.gz"))[, c("chunk", "chr", "st", "sp", "PPA_3")]
	colnames(chunkFile) = c("chunk", "chunk_CHR", "chunk_ST", "chunk_SP", "chunkPPA_3")

	## maximum snpPPA3 for each chunk (no chunks were filtered)
	merged = merge(chunkFile, max_snpPPA, by = "chunk")
	merged = merged[order(merged$chunkPPA_3, merged$snpPPA_3, decreasing = TRUE), ]
	output_noFilter = paste0(save_path, ".noFilter.txt.gz")
	fwrite(merged, file = output_noFilter, sep="\t", quote=FALSE)
	
	## filter chunk with regional_PPA3>0.9
	merged_Filter = merged[merged$chunkPPA_3>=chunkPPA3 & merged$snpPPA_3>=snpPPA3, ] 
	output_Filter = paste0(save_path, "_", chunkPPA3, "chunkPPA_", snpPPA3, "snpPPA.txt.gz")
	fwrite(merged_Filter, file = output_Filter, sep="\t", quote=FALSE)
	
	return(list(nofilter=output_noFilter, sig=output_Filter))
}

