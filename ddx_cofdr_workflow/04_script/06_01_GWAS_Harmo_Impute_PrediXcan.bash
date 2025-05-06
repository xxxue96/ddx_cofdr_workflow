# !/bin/bash

# gwas harmonization and imputation based on european 1000 Genomes hg38 reference panel https://github.com/hakyimlab/summary-gwas-imputation/wiki/GWAS-Harmonization-And-Imputation
# $1: phenotype
# $2: dir to save gwas
# $3: dir to save harmonized gwas
# export GWAS_TOOLS=/mnt/data/xue/Tool/summary-gwas-imputation/src
# export DATA=/mnt/data/xue/Tool/data_MetaXcan_hg38
source export_params.bash

pheno=$1
INPUT=$2
OUTPUT=$3

# step1: harmonize
## for gwas without FRQ: standardise original gwas with MungeSumstats -> harmo_gwas -> standardise again with MungeSumstats
harmo_gwas()
{
  echo "$pheno harmonization start"
  HARMO_FILE=$OUTPUT/01_harmonized_gwas/$pheno.txt
  
  if [ -f $HARMO_FILE.gz ]; then
    echo "$HARMO_FILE.gz exists."
  else
    python $GWAS_TOOLS/gwas_parsing.py \
    -gwas_file $INPUT/$pheno.txt.gz \
    -liftover $DATA/liftover/hg19ToHg38.over.chain.gz \
    -snp_reference_metadata $DATA/reference_panel_1000G/variant_metadata.txt.gz METADATA \
    -output_column_map SNP variant_id \
    -output_column_map CHR chromosome \
    -output_column_map BP position \
    -output_column_map A1 non_effect_allele \
    -output_column_map A2 effect_allele \
    -output_column_map FRQ frequency \
    -output_column_map BETA effect_size \
    -output_column_map SE standard_error \
    -output_column_map P pvalue \
    -output_column_map Z zscore \
    -output_column_map N sample_size \
    --chromosome_format \
    --enforce_numeric_columns \
    -output_order variant_id panel_variant_id chromosome position effect_allele non_effect_allele frequency pvalue zscore effect_size standard_error sample_size\
    -output $HARMO_FILE

    ## remove rows with NA
    sed -i '/NA/d' $HARMO_FILE && echo "$(wc -l $HARMO_FILE) snps after removing NA"
    gzip $HARMO_FILE	
  fi
  
  echo "$pheno harmonization done"
}

# step2: imputation
## --frequency_filter 0.01: only common variants included
impute_gwas()
{
  echo "$pheno imputation start"
  
  for i in {1..22}
  do
    for j in {0..9}
	do
	  echo "chr$i sub_batch$j imputation start"
	  
	  python $GWAS_TOOLS/gwas_summary_imputation.py \
	  -by_region_file $DATA/eur_ld.bed.gz \
	  -gwas_file $OUTPUT/01_harmonized_gwas/$pheno.txt.gz \
	  -parquet_genotype $DATA/reference_panel_1000G/chr$i.variants.parquet \
	  -parquet_genotype_metadata $DATA/reference_panel_1000G/variant_metadata.parquet \
	  -window 100000 \
	  -parsimony 7 \
	  -chromosome $i \
	  -regularization 0.1 \
	  -frequency_filter 0.01 \
	  -sub_batches 10 \
	  -sub_batch $j \
	  --standardise_dosages \
	  --cache_variants \
	  -output $OUTPUT/02_summary_imputation/"$pheno"_chr"$i"_sb"$j"_reg0.1_ff0.01_by_region.txt.gz
	  
	  echo "chr$i sub_batch$j imputation done"
	done
	
	echo "$pheno imputation done"
  done
}

# step3: imputation post processing
## gather all sun_batch analysis results
impute_post_process()
{
  echo "$pheno imputation postprocess start"
  
  python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
  -gwas_file $OUTPUT/01_harmonized_gwas/$pheno.txt.gz\
  -folder $OUTPUT/02_summary_imputation \
  -pattern $pheno.* \
  -parsimony 7 \
  -output $OUTPUT/03_processed_summary_imputation/"$pheno".imputed.txt.gz
  
  rm $OUTPUT/01_harmonized_gwas/$pheno.txt.gz
  rm $OUTPUT/02_summary_imputation/$pheno*
  echo "$pheno imputation postprocess done"
}
