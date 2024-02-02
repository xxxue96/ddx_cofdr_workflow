# ddx_cofdr_workflow
Differentiation and colocalization analysis at SNP- and gene-levels
## General description
1. Data preprocessing
2. Differentiation analysis
3. Colocalization analysis
## Working environment
python=3.7.12

R=4.2.1

Download scripts to your local working directory.
```
unzip ddx_cofdr_workflow.zip && cd ddx_cofdr_workflow
# Pleases modify paths of required files/tools or other paramaters in params.r and export_params.bash before any analysis.
```
## 1. Data preprocessing
Make subdirectories
```
mkdir -p 01_data/{00_Standardise_GWAS,01_LDSC_Results/{sumstats,h2,rg},02_tmvtnorm_Results,03_Gene_Result,04_FUMA_Results/{SNP2GENE,Cell_Type}}
```
Starndardise GWAS in R and save it to 01_data/00_Standardise_GWAS. 
```
source("01_data/01_Standardise_GWAS.r")
gwas = fread(raw_gwas_path) %>% mutate(FRQ = as.numeric(FRQ)) %>% select(SNP, CHR, BP, EFFECT_ALLELE, NON_EFFECT_ALLELE, FRQ, BETA, SE, P, Z, N, INFO)
gwas_formatted = Standardise_GWAS(gwas_sumstats_path=gwas, save_path=paste0("01_data/00_Standardise_GWAS/", trait, ".txt.gz"), ref_genome="GRCH37", convert_ref_genome=NULL, snp_ids_are_rs_ids = TRUE, INFO_filter = 0.3, FRQ_filter = 0.01, return_path = TRUE)
```
Run LDSC, and save rg and intercept with traits as rownames and colnames to 01_data/LDSC_Results. 
```
export ldsc=ldsc_tool_directory # download from https://github.com/bulik/ldsc
bash 04_script/02_01_Run_LDSC.bash t1 t2 # multiple traits are accepted
Rscript 04_script/02_02_Extract_LDSC.r t1 t2
```
Run mle.tmvnorm to estimate correlation coefficient due to sample overlap, and save estimates with traits as rownames and colnames to 01_data/02_tmvtnorm_Results. (For traits with invalidated LDSC estimates due to limited sample size)
```
Rscript 04_script/02_03_Fit_tmvtnorm.r t1 t2
```
Run MAGMA, TWAS, spTWAS, and PWAS for each trait in R, and save results to 01_data/03_Gene_Result.
```
rm(list=ls())
source("02_ddx/05_01_Run_DDx_Gene_MAGMA.r")
source("02_ddx/05_02_Run_DDx_Gene_PrediXcan.r")
source("02_ddx/05_03_Run_DDx_Gene_PWAS.r")
source("params.r") # please modify paramaters manually
FDRthres = 1 # keep all MAGMA/TWAS/spTWAS/PWAS-detected genes for further annotation
lapply(c(t1, t2), function(pheno){
	ddx_gwas = paste0("01_data/00_Standardise_GWAS/", pheno, ".txt.gz")
	twas_dir_ddx = paste0("01_data/03_Gene_Result/", pheno)
	magma_dir_ddx = paste0(twas_dir_ddx, "/01_MAGMA")
	pwas_dir_ddx = paste0(twas_dir_ddx, "/04_pwas")
	sapply(list(twas_dir_ddx, magma_dir_ddx, pwas_dir_ddx), dir.create, showWarnings = FALSE)

	Run_DDx_Gene_MAGMA(decor_gwas_paths=ddx_gwas, magma_dir=magma_dir_ddx, save_dir=twas_dir_ddx)
	Run_DDx_Gene_PrediXcan(traits=pheno, decor_gwas_paths=ddx_gwas, metaxcan_dir=twas_dir_ddx, save_dir=twas_dir_ddx)
	Run_DDx_Gene_PWAS(traits=pheno, decor_gwas_paths=ddx_gwas, pwas_dir=pwas_dir_ddx, save_dir=twas_dir_ddx)

	# remove temporary files
	file.remove(paste0(twas_dir_ddx,"/01_harmonized_gwas/",pheno,".txt.gz"))
})
```
Analyze gwas via FUMA https://fuma.ctglab.nl "SNP2GENE" and "Cell Type" module. Save corresponding zipped results files to 01_data/04_FUMA_Results/{SNP2GENE,Cell_Type} and rename results files with trait names

## 2. Differeniation analysis
Run SNP- and gene-level differentiation analysis
```
source("02_ddx/DDx_Workflow.r")
source("02_ddx/DDx_Workflow_params.r")
# snp-level
ddx_wrap_up(t1, t2)

# Analyze ddx-derived gwas via FUMA, and save results to 02_ddx/04_FUMA_Results/{SNP2GENE,Cell_Type,Image}

# gene-level
ddx_wrap_up_gene(t1, t2)
```
Format results to readable excel files
```
source("02_ddx/06_01_DDx_Format_Excel.r")
source("02_ddx/06_02_DDx_Format_Image.r")
source("02_ddx/06_03_DDx_Excel_Summary.r")
source("02_ddx/06_03_DDx_Excel_Summary_more_info.r")
source("02_ddx/06_04_DDx_Count_Res.r")
# format differentiation results to excel
DDx_Format_Excel(t1, t2, save_dir="02_ddx/13_Excel_Summary/Tables")

# summarize all images from FUMA SNP2GENE and Cell Type module (multiple t1 and t2 are accepted)
DDx_Format_Image(t1, t2, save_path="02_ddx/13_Excel_Summary/Figures/all_trait_comparison_image.xlsx")

# summarize all excel files in excel_dir (multiple excel_dir are accepted)
DDx_Excel_Summary(excel_dir="02_ddx/13_Excel_Summary/Tables", save_path="02_ddx/13_Excel_Summary/Tables//all_trait_comparison_summary.xlsx")
DDx_Excel_Summary_more_info(excel_dir="02_ddx/13_Excel_Summary/Tables", save_path="02_ddx/13_Excel_Summary/Tables/all_trait_comparison_summary_more_info.xlsx")

# summarize number of significant results 
DDx_Count_Res(excel_dir="02_ddx/13_Excel_Summary/Tables", save_path="02_ddx/13_Excel_Summary/all_trait_comparison_DDx_Count.txt")
```
## 3. Colocalization analysis
Run SNP- and gene-level colocalization analysis
```
source("03_cofdr/Cofdr_Workflow.r")
source("03_cofdr/Cofdr_Workflow_params.r")
# snp-level
cofdr_wrap_up(t1, t2)

# Analyze cofdr-derived gwas via FUMA, and save results to 03_cofdr/05_02_Cofdr_FUMA_Result/{SNP2GENE,GENE2FUNC}

# gene-level
cofdr_wrap_up_gene(t1, t2)
```
Format results to readable excel files
```
source("03_cofdr/06_01_Cofdr_Format_Excel.r")
source("03_cofdr/06_03_Cofdr_Excel_Summary.r")
source("03_cofdr/06_03_Cofdr_Excel_Summary_more_info.r")
source("03_cofdr/06_04_Cofdr_Count_Res.r")
# format colocalization results to excel
Cofdr_Format_Excel(t1, t2, save_dir="03_cofdr/10_Excel_Summary/Tables")

# summarize all images from FUMA GENE2FUNC module (multiple t1 and t2 are accepted)
Cofdr_Format_Image(t1, t2, save_path = "03_cofdr/10_Excel_Summary/Figures/all_trait_comparison_image.xlsx")

# summarize all excel files in excel_dir (multiple excel_dir are accepted)
Cofdr_Excel_Summary(excel_dir="03_cofdr/10_Excel_Summary/Tables", save_path="03_cofdr/10_Excel_Summary/all_trait_comparison_summary.xlsx")
Cofdr_Excel_Summary_more_info(excel_dir="03_cofdr/10_Excel_Summary/Tables", save_path="03_cofdr/10_Excel_Summary/all_trait_comparison_summary_more_info.xlsx")

# summarize number of significant results
Cofdr_Count_Res(excel_dir="03_cofdr/10_Excel_Summary/Tables/", save_path="03_cofdr/10_Excel_Summary/all_trait_comparison_Cofdr_Count.txt")
```






