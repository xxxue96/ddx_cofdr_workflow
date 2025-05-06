# Run LDSC
# conda activate py27
# export ldsc=/mnt/data/xue/Tool/ldsc
source export_params.bash

traits=$@
trait_length=$#

echo "Run LDSC among traits: $traits"
echo "Number of traits: $trait_length"

# select snps overlap with hapmap3 snps for ldsc regresion elements, abd w_hm3.snplist is well-imputed snps in most studies 
# a1/a2 alleles in munge_sumstat.py refers to effect/non-effect alleles https://github.com/bulik/ldsc/issues/133, which referes to A2/A1 in standardized gwas from library(MungeSumstats)
# ignore columns that may contain NA
for trait in $traits
do
  python $ldsc/munge_sumstats.py --sumstats 01_data/00_Standardise_GWAS/"$trait".txt.gz --out 01_data/01_LDSC_Results/sumstats/"$trait" --merge-alleles $ldsc/ld_ref/hapmap3/w_hm3.snplist --chunksize 500000 \
  --N-col N \
  --a1 A2 --a2 A1 \
  --signed-sumstats Z,0 \
  --ignore BETA,SE,FRQ 
done

# cal. observed scale h2 when continous traits
for trait in $traits
do
  echo "cal. observed scale h2"
  python $ldsc/ldsc.py --h2 01_data/01_LDSC_Results/sumstats/"$trait".sumstats.gz --ref-ld-chr $ldsc/ld_ref/eur_w_ld_chr/ --w-ld-chr $ldsc/ld_ref/eur_w_ld_chr/ --out 01_data/01_LDSC_Results/h2/"$trait"_h2
done

# cal. rg
for (( i=1; i<=$trait_length; i++ ))
do
  for ((j=i+1; j<=$trait_length; j++))
  do
    if [ ! -f "01_data/01_LDSC_Results/rg/${!i}_${!j}.log" ]; then
	  python $ldsc/ldsc.py --rg 01_data/01_LDSC_Results/sumstats/${!i}.sumstats.gz,01_data/01_LDSC_Results/sumstats/${!j}.sumstats.gz --ref-ld-chr $ldsc/ld_ref/eur_w_ld_chr/ --w-ld-chr $ldsc/ld_ref/eur_w_ld_chr/ --out 01_data/01_LDSC_Results/rg/${!i}_${!j}
	fi
  done
done


# rg is nan when h2 is negative, then use phenotypic correlation or fit tmvtnorm to estimate sample overlap correlation instead of ldsc intercept

# cal. liability h2 when binary traits (assuming poppulation prevenlence = sample prevenlence = N_CAS/N)
for trait in $traits
do
  echo "cal. liability h2"
  export t1_prev=$(awk -F'\t' -v abbre=$trait '$1 == abbre {print $7}' $sample_size_src)
  python $ldsc/ldsc.py --h2 01_data/01_LDSC_Results/sumstats/"$trait".sumstats.gz --ref-ld-chr $ldsc/ld_ref/eur_w_ld_chr/ --w-ld-chr $ldsc/ld_ref/eur_w_ld_chr/ --out 01_data/01_LDSC_Results/h2/"$trait"_h2lia --samp-prev $t1_prev --pop-prev $t1_prev
done

# cal. rg
# for (( i=1; i<=$trait_length; i++ ))
# do
  # for ((j=i+1; j<=$trait_length; j++))
  # do
    # if [ ! -f "01_data/01_LDSC_Results/rg/${!i}_${!j}_liability.log" ]; then
	  # export t1_prev=$(awk -F'\t' -v abbre=${!i} '$1 == abbre {print $7}' $sample_size_src)
      # export t2_prev=$(awk -F'\t' -v abbre=${!j} '$1 == abbre {print $7}' $sample_size_src)
	  # python $ldsc/ldsc.py --rg 01_data/01_LDSC_Results/sumstats/${!i}.sumstats.gz,01_data/01_LDSC_Results/sumstats/${!j}.sumstats.gz --ref-ld-chr $ldsc/ld_ref/eur_w_ld_chr/ --w-ld-chr $ldsc/ld_ref/eur_w_ld_chr/ --out 01_data/01_LDSC_Results/rg/${!i}_${!j}_liability --samp-prev $t1_prev,$t2_prev --pop-prev $t1_prev,$t2_prev
	# fi
  # done
# done
