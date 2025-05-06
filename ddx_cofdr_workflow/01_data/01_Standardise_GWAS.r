# https://github.com/neurogenomics/MungeSumstats/blob/master/vignettes/MungeSumstats.Rmd

# Standardise gwas summary statistics as the format below. 
# SNP : SNP ID (rs IDs)
# CHR : Chromosome number
# BP : Base pair positions
# A1 : reference allele and non-effect allele
# A2 : alternative allele and effect allele
# Z : Z-score (by default the Z-score is assumed to be calculated off the effect sample_size not the P-value and so will be flipped if necessary. This can be changed by a user. it can be imputed from P value with "compute_z = TRUE")
# BETA : Effect sample_size estimate relative to the alternative allele
# P : Unadjusted p-value for SNP
# SE : The standard error
# N : Sample sample_size
# INFO: The imputation information score
# FRQ: Frequency of A2 i.e. EAF 
#	   Converntionally, it is also MAF since A2 usually is also minor allele; Sometimes, A2 can be major allele and FRQ will be renamed as "MAJOR_ALLELE_FRQ" if set "frq_is_maf = FALSE"

# data("sumstatsColHeaders") to check recoginized colmun headers
# Input can be data.table/data.frame or file paths of VCF, txt, tsv, csv file types or .gz/.bgz versions of these file types
# Output can be data.table/data.frame or file paths of reformatted gwas

# Reference genome can be specified by user or infer from raw gwas (ref_genome=NULL), download link:
# BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
# BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
# BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
# BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

# conda activate R42

library(MungeSumstats)
library(data.table)
setDTthreads(40)

Standardise_GWAS <- function(gwas_sumstats_path, save_path, ref_genome = "GRCH37", convert_ref_genome = NULL, snp_ids_are_rs_ids = TRUE, INFO_filter = 0.3, FRQ_filter = 0, return_path = TRUE)
{
	data("sumstatsColHeaders")
	gwas = format_sumstats( 
							path = gwas_sumstats_path,
							ref_genome = ref_genome, # GRCH37, GRCH38, or NULL if inferring from gwas_sumstats_path
							convert_ref_genome = convert_ref_genome, # GRCH37, GRCH38, or NULL if no ref_genome conversion
							convert_small_p = TRUE,
							convert_large_p = TRUE,
							convert_neg_p = TRUE,
							compute_z = FALSE,
							force_new_z = FALSE,
							compute_n = 0L,
							convert_n_int = TRUE,
							impute_beta = FALSE,
							impute_se = FALSE,
							analysis_trait = NULL,
							INFO_filter = INFO_filter, # filter snps with imputation quality if there is an INFO column
							FRQ_filter = FRQ_filter, # filter snps with MAF
							pos_se = TRUE,
							effect_columns_nonzero = FALSE,
							N_std = 5,
							N_dropNA = TRUE,
							rmv_chr = c("X", "Y", "MT"), 
							rmv_chrPrefix = TRUE,
							on_ref_genome = TRUE,
							strand_ambig_filter = FALSE, 
							allele_flip_check = TRUE,
							allele_flip_drop = TRUE,
							allele_flip_z = TRUE,
							allele_flip_frq = TRUE,
							bi_allelic_filter = TRUE,
							snp_ids_are_rs_ids = snp_ids_are_rs_ids, # FALSE if imputing rsid from chr:bp
							remove_multi_rs_snp = FALSE,
							frq_is_maf = TRUE,
							indels = TRUE,
							dbSNP = 144, # version of dbSNP to be used for imputation (144vs155 https://www.biostars.org/p/9534249/)
							sort_coordinates = TRUE,
							nThread = 40,
							save_path = tempfile(fileext = ".tsv.gz"),
							write_vcf = FALSE,
							tabix_index = FALSE,
							return_data = TRUE,
							return_format = "data.table",
							ldsc_format = FALSE,
							save_format = NULL,
							log_folder_ind = FALSE,
							log_mungesumstats_msgs = FALSE,
							imputation_ind = FALSE,
							log_folder = tempdir(),
							force_new = FALSE,
							mapping_file = sumstatsColHeaders
						  )

	gwas = gwas[!duplicated(gwas$SNP), c("SNP", "CHR", "BP", "A1", "A2", "FRQ", "BETA", "SE", "P", "Z", "N")]	

	if(return_path == TRUE){
		fwrite(gwas, file = save_path, sep = "\t", na = "NA", quote = FALSE, scipen = 999)
		return(save_path)
	} else{
		return(gwas)
	}
}


