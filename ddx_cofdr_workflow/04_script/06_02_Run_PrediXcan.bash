#!/bin/bash

# integrating GWAS and QTL studies with methods predixcan and multixcan, i.e. map snps to genes based on qtl information
# https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS

# ----------------------------------------------------------------------------------------------------------------------#																		#
# 	GTEX_v8_MASHR_models: whether and how snps affect gene expression and splicing levels measured in GTEX_v8 database; #
#						  -> effect sizes computed with MASHR and fine-mapping analysis done by DAP-G 					#
#						  -> regression model with gene expression as outcome and snp dosage as independent variable	#
#	Harmonization and imputation: not all fine-mapped snps in GTEX_v8_MASHR_models also appeared in our gwas; 			#
#								  -> harmonize and impute missing ine-mapped variants' summary statistics from the GWAS	#
#									 both based on 1kg reference panel;													#
#									 GTEx data as genotype reference panel would be better, but the data is non-shared	#														#																																		#
# ----------------------------------------------------------------------------------------------------------------------#

# $1: phenotype
# $2: dir to save gwas
# $3: dir to save harmonized gwas
# $4: dir to save metaxcan results
source 04_script/06_01_GWAS_Harmo_Impute_PrediXcan.bash

pheno=$1
INPUT=$2
OUTPUT=$3
OUTPUT_Res=$4

# step1: harmonization to hg38 1kg reference panel
if [[ "$run_harmo" == 1 ]]; then harmo_gwas; fi

# step2: imputation to hg38 1kg reference panel
## ignore it due to time consuming, but there are still 81 % of model's snps used so it's ok
if [[ "$run_impute" == 1 ]]; then impute_gwas; impute_post_process; fi

# step3: run spredixcan for 49 GTEx tissues
## $1: eqtl/sqtl
## $2: tissue name
run_spredixcan()
{
  mkdir -p $OUTPUT_Res/02_spredixcan/$1/$pheno
  if [[ "$run_impute" == 1 ]]; then GWAS_FILE=$OUTPUT/03_processed_summary_imputation/"$pheno"_impute.txt.gz; else GWAS_FILE=$OUTPUT/01_harmonized_gwas/$pheno.txt.gz; fi

  python $METAXCAN/SPrediXcan.py \
  --gwas_file $GWAS_FILE \
  --snp_column panel_variant_id \
  --effect_allele_column effect_allele \
  --non_effect_allele_column non_effect_allele \
  --zscore_column zscore \
  --model_db_path $DATA/models/$1/mashr/mashr_$2.db \
  --covariance $DATA/models/$1/mashr/mashr_$2.txt.gz \
  --keep_non_rsid \
  --additional_output \
  --model_db_snp_key varID \
  --throw \
  --output_file $OUTPUT_Res/02_spredixcan/$1/$pheno/"$pheno"__PM__$2.csv
}

# step4: run smultixcan
## $1: eqtl/sqtl
run_smultixcan()
{
  mkdir -p $OUTPUT_Res/03_smultixcan/$1
  if [[ "$run_impute" == 1 ]]; then GWAS_FILE=$OUTPUT/03_processed_summary_imputation/"$pheno"_impute.txt.gz; else GWAS_FILE=$OUTPUT/01_harmonized_gwas/$pheno.txt.gz; fi
  if [[ "$1" == "eqtl" ]]; then SNP_COVAR=$DATA/models/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz; else SNP_COVAR=$DATA/models/gtex_v8_splicing_mashr_snp_smultixcan_covariance.txt.gz; fi

  python $METAXCAN/SMulTiXcan.py \
  --models_folder $DATA/models/$1/mashr \
  --models_name_pattern "mashr_(.*).db" \
  --models_name_filter "mashr_(.*).db" \
  --snp_covariance $SNP_COVAR \
  --metaxcan_folder $OUTPUT_Res/02_spredixcan/$1/$pheno \
  --metaxcan_filter "$pheno__PM__(.*).csv" \
  --metaxcan_file_name_parse_pattern "(.*)__PM__(.*).csv" \
  --gwas_file $GWAS_FILE \
  --snp_column panel_variant_id \
  --effect_allele_column effect_allele \
  --non_effect_allele_column non_effect_allele \
  --zscore_column zscore \
  --keep_non_rsid \
  --model_db_snp_key varID \
  --cutoff_condition_number 30 \
  --verbosity 7 \
  --throw \
  --output $OUTPUT_Res/03_smultixcan/$1/"$pheno"_smultixcan.txt
}

run_smultixcan()
{
  mkdir -p $OUTPUT_Res/03_smultixcan/$1
  if [[ "$run_impute" == 1 ]]; then GWAS_FILE=$OUTPUT/03_processed_summary_imputation/"$pheno"_impute.txt.gz; else GWAS_FILE=$OUTPUT/01_harmonized_gwas/$pheno.txt.gz; fi
  if [[ "$1" == "eqtl" ]]; then SNP_COVAR=$DATA/models/gtex_v8_expression_mashr_snp_smultixcan_covariance.txt.gz; else SNP_COVAR=$DATA/models/gtex_v8_splicing_mashr_snp_smultixcan_covariance.txt.gz; fi

  python $METAXCAN/SMulTiXcan.py \
  --models_folder $DATA/models/$1/mashr \
  --models_name_pattern "mashr_(.*).db" \
  --models_name_filter $(printf "mashr_(%s).db" $(tr '\n' '|' < $model_src | sed 's/|$//')) \
  --snp_covariance $SNP_COVAR \
  --metaxcan_folder $OUTPUT_Res/02_spredixcan/$1/$pheno \
  --metaxcan_filter "$pheno__PM__(.*).csv" \
  --metaxcan_file_name_parse_pattern "(.*)__PM__(.*).csv" \
  --gwas_file $GWAS_FILE \
  --snp_column panel_variant_id \
  --effect_allele_column effect_allele \
  --non_effect_allele_column non_effect_allele \
  --zscore_column zscore \
  --keep_non_rsid \
  --model_db_snp_key varID \
  --cutoff_condition_number 30 \
  --verbosity 7 \
  --throw \
  --output $OUTPUT_Res/03_smultixcan/$1/"$pheno"_smultixcan.txt
}

# run s-predixcan for each tissue
for ts in $(cat $model_src)
do
  run_spredixcan eqtl $ts
  run_spredixcan sqtl $ts
done

# run s-multixcan across multiple tissues 
run_smultixcan eqtl 
run_smultixcan sqtl 
