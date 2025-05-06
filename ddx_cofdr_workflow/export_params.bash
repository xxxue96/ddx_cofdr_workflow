# please modify paths of required files/tools or other paramaters before any analysis 

# 1. tools to be installed 
# LDSC
export ldsc=/mnt/data/xue/Tool/ldsc # https://github.com/bulik/ldsc

# MetaXcan https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS
export GWAS_TOOLS=/mnt/data/xue/Tool/summary-gwas-imputation/src # https://github.com/hakyimlab/summary-gwas-imputation/tree/master
export DATA=/mnt/data/xue/Tool/data_MetaXcan_hg38 # https://zenodo.org/records/3657902#.Xj2Zh-F7m90
export METAXCAN=/mnt/data/xue/Tool/MetaXcan/software # https://github.com/hakyimlab/MetaXcan

# PWAS
export PWAS=/mnt/data/xue/Tool/PWAS # http://nilanjanchatterjeelab.org/pwas/

# 2. files to be downloaded

# 3. other files to be prepared
# sample size with "N_CAS" as 2nd column
export sample_size_src=01_data/sample_size.txt

# PrediXcan models with tissues as 1st column
export model_src=05_auxiliary/predixcan_masher_models.txt

# whether run harmonization in 04_script/06_02_Run_PrediXcan.bash
export run_harmo=1

# whether run imputation in 04_script/06_02_Run_PrediXcan.bash (time consuming, here we ignore this step)
export run_impute=0



