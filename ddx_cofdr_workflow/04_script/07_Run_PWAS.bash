# download FUSION predictive models for plasma proteins
# model1. https://interval-pwas.s3.us-west-1.amazonaws.com/INTERVAL.FUSION.tar.gz PMID: 34315903; N=3301EUR
# model2. http://nilanjanchatterjeelab.org/pwas/ PMID: 35501419; N=7213EA / N=1871AA (hg19 recommended)

# export PWAS=/mnt/data/xue/Tool/PWAS
source export_params.bash

# $1: phenotype
# $2: gwas summary statistics, A1 as effect allele
# $3: directory to save analyzed PWAS results
pheno=$1
INPUT=$2
OUTPUT=$3

CHR=22
until [ $CHR -lt 1 ]
do  
  if [ -f "$OUTPUT/"$pheno"_chr${CHR}.out" ]; then
    echo "$OUTPUT/"$pheno"_chr${CHR}.out exists."
  else
    Rscript $PWAS/scripts/PWAS.assoc_test.R \
    --sumstats $INPUT \
    --weights $PWAS/PWAS_EA/Plasma_Protein_EA_hg19.pos \
    --weights_dir $PWAS/PWAS_EA/Plasma_Protein_weights_EA/ \
    --ref_ld_chr $PWAS/LDref/EUR/chr \
    --force_model enet \
    --chr ${CHR} \
    --out $OUTPUT/"$pheno"_chr${CHR}.out  
  fi 
  let CHR-=1
done
